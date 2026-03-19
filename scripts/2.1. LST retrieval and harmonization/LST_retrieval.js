var CITY_KEY = 'IASI'; 
var CITIES = {
  'BUCURESTI':    {name:'București',    lon:26.1025, lat:44.4268, zoom:10},
  'CLUJ_NAPOCA': {name:'Cluj-Napoca', lon:23.5899, lat:46.7712, zoom:11},
  'IASI':         {name:'Iași',         lon:27.5879, lat:47.1585, zoom:11},
  'CONSTANTA':    {name:'Constanța',    lon:28.6348, lat:44.1598, zoom:11},
  'TIMISOARA':    {name:'Timișoara',    lon:21.2287, lat:45.7489, zoom:11},
  'BRASOV':       {name:'Brașov',       lon:25.6012, lat:45.6579, zoom:11},
  'CRAIOVA':      {name:'Craiova',      lon:23.7949, lat:44.3302, zoom:11}
};

var city = CITIES[CITY_KEY];
if (!city) throw new Error('Unknown CITY_KEY: ' + CITY_KEY);

var WORK_RADIUS_KM = 35;
var AOI = ee.Geometry.Point([city.lon, city.lat]).buffer(WORK_RADIUS_KM * 1000);

Map.setCenter(city.lon, city.lat, city.zoom);
Map.addLayer(AOI, {}, 'AOI');

var YEARS = [2021, 2022, 2023, 2024, 2025];
var MONTH_START = 6; 
var MONTH_END   = 8;

var MIN_INTERSECTION_FRAC = 0; 
var MAX_CLOUD_FRAC_IN_FOOTPRINT = 0.10; 
var STATS_SCALE = 120;

var DRIVE_FOLDER_RASTER = 'GEE_LST_yearly_HARMONIZED';
var DRIVE_FOLDER_CSV    = 'GEE_LST_dailyStats_HARMONIZED_CSV';
var OUT_CRS      = 'EPSG:3035';
var OUT_SCALE_M  = 30;
var MAX_PIXELS   = 1e13;
var FILL_NODATA_FINAL = false;
var NODATA_VALUE = -9999;

function prepL2(img) {
  var qa = img.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1<<0).eq(0).and(qa.bitwiseAnd(1<<3).eq(0));
  
  var lst = img.select('ST_B10')
    .multiply(0.00341802).add(149.0).subtract(273.15) 
    .rename('LST_C');
    
  return lst.updateMask(mask)
    .addBands(qa) 
    .copyProperties(img, ['system:time_start', 'SPACECRAFT_ID']);
}

function prepTOA(img) {
  var qa = img.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1<<0).eq(0).and(qa.bitwiseAnd(1<<3).eq(0));
  
  var toa = img.select('B10')
    .subtract(273.15) 
    .rename('TOA');
    
  return toa.updateMask(mask)
    .copyProperties(img, ['system:time_start']);
}

function buildHarmonizedDailyMosaics(year) {
  var start = ee.Date.fromYMD(year, MONTH_START, 1);
  var end   = ee.Date.fromYMD(year, MONTH_END + 1, 1); 

  var l8_l2 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2');
  var l9_l2 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2');
  var col_L2 = l8_l2.merge(l9_l2)
    .filterBounds(AOI)
    .filterDate(start, end)
    .map(prepL2);

  var l8_toa = ee.ImageCollection('LANDSAT/LC08/C02/T1_TOA');
  var l9_toa = ee.ImageCollection('LANDSAT/LC09/C02/T1_TOA');
  var col_TOA = l8_toa.merge(l9_toa)
    .filterBounds(AOI)
    .filterDate(start, end)
    .map(prepTOA);

  var dates = ee.List(col_L2.map(function(img){
    return ee.Feature(null, {'dateKey': ee.Date(img.get('system:time_start')).format('YYYYMMdd')});
  }).aggregate_array('dateKey')).distinct().sort();

  var dailyMosaics = dates.map(function(dateStr) {
    dateStr = ee.String(dateStr);
    var dayStart = ee.Date.parse('YYYYMMdd', dateStr);
    var dayEnd = dayStart.advance(1, 'day');

    var dailyL2 = col_L2.filterDate(dayStart, dayEnd);
    var dailyTOA = col_TOA.filterDate(dayStart, dayEnd);

    var filter = ee.Filter.equals({leftField: 'system:time_start', rightField: 'system:time_start'});
    var joined = ee.Join.saveFirst('toa_match').apply(dailyL2, dailyTOA, filter);

    var processedScenes = ee.ImageCollection(joined.map(function(imgL2) {
      imgL2 = ee.Image(imgL2);
      var imgTOA = ee.Image(imgL2.get('toa_match'));
      
      return ee.Algorithms.If(imgTOA, 
        (function(){
          var regInput = imgTOA.addBands(imgL2).select(['TOA', 'LST_C']);
          var fit = regInput.reduceRegion({
            reducer: ee.Reducer.linearFit(),
            geometry: AOI,
            scale: 120, 
            bestEffort: true,
            maxPixels: 1e9
          });
          
          var scale = ee.Number(fit.get('scale'));
          var offset = ee.Number(fit.get('offset'));
          
          scale = ee.Algorithms.If(scale, scale, 1);
          offset = ee.Algorithms.If(offset, offset, 0);
          
          var scaleImg = ee.Image.constant(scale);
          var offsetImg = ee.Image.constant(offset);
          
          var predLST = imgTOA.select('TOA')
            .multiply(scaleImg)
            .add(offsetImg)
            .rename('LST_C');
          
          var patched = imgL2.select('LST_C').unmask(predLST);
          
          return patched
            .addBands(imgL2.select('QA_PIXEL')) 
            .copyProperties(imgL2, ['system:time_start']);
        })(),
        imgL2
      );
    }));

    var mosaicLST = processedScenes.select('LST_C').mean(); 
    var mosaicQA  = processedScenes.select('QA_PIXEL').first(); 

    var obsAvail = mosaicLST.mask().rename('obsAvail'); 
    var obs01 = obsAvail.unmask(0);

    var intersectionFrac = ee.Number(
      obs01.reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: AOI,
        scale: STATS_SCALE,
        bestEffort: true,
        maxPixels: MAX_PIXELS
      }).get('obsAvail')
    );

    var cloudBit = mosaicQA.bitwiseAnd(1<<3).neq(0);
    var clearBit = cloudBit.not();
    
    var clearFrac = ee.Number(
      clearBit.updateMask(obsAvail).reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: AOI,
        scale: STATS_SCALE,
        bestEffort: true,
        maxPixels: MAX_PIXELS
      }).get('QA_PIXEL')
    );
    
    clearFrac = ee.Number(ee.Algorithms.If(clearFrac, clearFrac, 0));
    var cloudFrac = ee.Number(1).subtract(clearFrac);

    return mosaicLST
      .clip(AOI)
      .set({
        city_key: CITY_KEY,
        year: year,
        date: dayStart.format('YYYY-MM-dd'),
        dateKey: dateStr,
        aoi_intersection_frac: intersectionFrac,
        aoi_cloud_frac_in_fp: cloudFrac
      });
  });

  return ee.ImageCollection.fromImages(dailyMosaics);
}


function buildSceneValidationFC(year) {
  var start = ee.Date.fromYMD(year, MONTH_START, 1);
  var end   = ee.Date.fromYMD(year, MONTH_END + 1, 1);

  var col_L2 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
    .filterBounds(AOI)
    .filterDate(start, end)
    .map(prepL2);

  var col_TOA = ee.ImageCollection('LANDSAT/LC08/C02/T1_TOA')
    .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_TOA'))
    .filterBounds(AOI)
    .filterDate(start, end)
    .map(prepTOA);

  var filter = ee.Filter.equals({
    leftField: 'system:time_start',
    rightField: 'system:time_start'
  });

  var joined = ee.Join.saveFirst('toa_match').apply(col_L2, col_TOA, filter);

  var fc = ee.FeatureCollection(joined.map(function(imgL2) {
    imgL2 = ee.Image(imgL2);
    var imgTOA = ee.Image(imgL2.get('toa_match'));

    return ee.Algorithms.If(imgTOA, (function() {
      var obs = imgL2.select('LST_C');
      var toa = imgTOA.select('TOA'); 

      var regInput = toa.addBands(obs).select(['TOA', 'LST_C']);
      var fit = regInput.reduceRegion({
        reducer: ee.Reducer.linearFit(),
        geometry: AOI,
        scale: 120,
        bestEffort: true,
        maxPixels: 1e9
      });

      var scale = ee.Number(fit.get('scale'));
      var offset = ee.Number(fit.get('offset'));

      scale  = ee.Number(ee.Algorithms.If(scale,  scale,  1));
      offset = ee.Number(ee.Algorithms.If(offset, offset, 0));

      var pred = toa.multiply(ee.Image.constant(scale))
                    .add(ee.Image.constant(offset))
                    .rename('pred');

      var resid = pred.subtract(obs).rename('resid');
      var sq    = resid.pow(2).rename('sq');
      var ab    = resid.abs().rename('abs');

      var obsMask  = obs.mask();
      var predMask = pred.mask();
      var fill01u = obsMask.not().and(predMask).rename('fill01').unmask(0);

      var stack = resid.addBands(sq).addBands(ab).addBands(fill01u);

      var stats = stack.reduceRegion({
        reducer: ee.Reducer.mean()
          .combine({reducer2: ee.Reducer.sum(),   sharedInputs: true})
          .combine({reducer2: ee.Reducer.count(), sharedInputs: true}),
        geometry: AOI,
        scale: STATS_SCALE,
        bestEffort: true,
        maxPixels: MAX_PIXELS
      });
      
      var biasObj = stats.get('resid_mean');
      var mseObj  = stats.get('sq_mean');
      var maeObj  = stats.get('abs_mean');
      var nObj    = stats.get('resid_count');

      var nFillObj   = stats.get('fill01_sum'); 
      var fillFracObj = stats.get('fill01_mean'); 
      
      var n = ee.Number(ee.Algorithms.If(nObj, nObj, 0));
      var bias = ee.Number(ee.Algorithms.If(biasObj, biasObj, -9999));
      var mae  = ee.Number(ee.Algorithms.If(maeObj,  maeObj,  -9999));
      
      var rmse = ee.Number(
        ee.Algorithms.If(
          n.gt(0),
          ee.Algorithms.If(mseObj, ee.Number(mseObj).sqrt(), -9999),
          -9999
        )
      );

      var n_fill = ee.Number(ee.Algorithms.If(nFillObj, nFillObj, 0));
      var fill_frac = ee.Number(ee.Algorithms.If(fillFracObj, fillFracObj, 0));

      var corrDict = pred.addBands(obs).reduceRegion({
        reducer: ee.Reducer.pearsonsCorrelation(),
        geometry: AOI,
        scale: STATS_SCALE,
        bestEffort: true,
        maxPixels: MAX_PIXELS
      });
      
      var rObj = corrDict.get('correlation');
      var r = ee.Number(ee.Algorithms.If(rObj, rObj, -9999));

      var t = ee.Date(imgL2.get('system:time_start'));

      return ee.Feature(null, {
        city_key: CITY_KEY,
        year: year,
        datetime: t.format('YYYY-MM-dd HH:mm'),
        date: t.format('YYYY-MM-dd'),
        system_time_start: imgL2.get('system:time_start'),
        spacecraft: imgL2.get('SPACECRAFT_ID'),
        scene_index: imgL2.get('system:index'),

        scale: scale,
        offset: offset,

        n_valid: n,
        bias_C: bias,
        mae_C: mae,
        rmse_C: rmse,
        pearson_r: r,

        n_fill: n_fill,
        fill_frac: fill_frac
      });
    })(), ee.Feature(null, {
      city_key: CITY_KEY,
      year: year,
      system_time_start: imgL2.get('system:time_start'),
      scene_index: imgL2.get('system:index'),
      note: 'No TOA match',
      n_fill: 0,
      fill_frac: 0
    }));
  }));

  return fc;
}


var multiYearCol = ee.ImageCollection([]);

YEARS.forEach(function(y) {
  print('Processing Year: ' + y + ' with Harmonization...');
  
  var dailyAll = buildHarmonizedDailyMosaics(y);

  var dailyKept = dailyAll
    .filter(ee.Filter.gte('aoi_intersection_frac', MIN_INTERSECTION_FRAC))
    .filter(ee.Filter.lte('aoi_cloud_frac_in_fp', MAX_CLOUD_FRAC_IN_FOOTPRINT));

  multiYearCol = multiYearCol.merge(dailyKept);

  print('Year', y, 'Days Total:', dailyAll.size(), 'Days Kept:', dailyKept.size());

  var annual = dailyKept.median().rename('LST_C');
  
  if (FILL_NODATA_FINAL) annual = annual.unmask(NODATA_VALUE);
  
  Map.addLayer(annual, {min: 15, max: 45, palette: ['blue', 'yellow', 'red']}, 'LST_Harmonized_' + y, false);

  var prefix = 'RO_' + CITY_KEY + '_LST_HARMONIZED_JJA_' + y +
    '_minInt' + Math.round(MIN_INTERSECTION_FRAC * 100) + 'pct' +
    '_maxCloudFP' + Math.round(MAX_CLOUD_FRAC_IN_FOOTPRINT * 100) + 'pct' +
    '_R' + WORK_RADIUS_KM + 'km';

  Export.image.toDrive({
    image: annual,
    description: prefix,
    folder: DRIVE_FOLDER_RASTER,
    fileNamePrefix: prefix,
    region: AOI,
    crs: OUT_CRS,
    scale: OUT_SCALE_M,
    maxPixels: MAX_PIXELS
  });

  var sceneValFC = buildSceneValidationFC(y);

  Export.table.toDrive({
    collection: sceneValFC,
    description: 'SCENE_VALIDATION_' + CITY_KEY + '_JJA_' + y,
    folder: DRIVE_FOLDER_CSV,
    fileNamePrefix: 'SCENE_VALIDATION_' + CITY_KEY + '_JJA_' + y,
    fileFormat: 'CSV'
  });

  var dailyStatsFC = ee.FeatureCollection(dailyAll.map(function(img){
      var inter = ee.Number(img.get('aoi_intersection_frac'));
      var cloud = ee.Number(img.get('aoi_cloud_frac_in_fp'));
      return ee.Feature(null, {
        city_key: CITY_KEY,
        year: y,
        date: img.get('date'),
        aoi_intersection_frac: inter,
        aoi_cloud_frac_in_fp: cloud,
        kept: inter.gte(MIN_INTERSECTION_FRAC).and(cloud.lte(MAX_CLOUD_FRAC_IN_FOOTPRINT))
      });
  }));
  
  Export.table.toDrive({
    collection: dailyStatsFC,
    description: 'DAILY_STATS_HARMONIZED_' + CITY_KEY + '_JJA_' + y,
    folder: DRIVE_FOLDER_CSV,
    fileNamePrefix: 'DAILY_STATS_HARMONIZED_' + CITY_KEY + '_JJA_' + y,
    fileFormat: 'CSV'
  });
});


print('Calculating Multi-Year Median from ' + multiYearCol.size() + ' days...');

var multiYearMedian = multiYearCol.median().rename('LST_C');

if (FILL_NODATA_FINAL) multiYearMedian = multiYearMedian.unmask(NODATA_VALUE);

var multiYearPrefix = 'RO_' + CITY_KEY + '_LST_HARMONIZED_JJA_MULTIYEAR_' + YEARS[0] + '-' + YEARS[YEARS.length-1] +
    '_minInt' + Math.round(MIN_INTERSECTION_FRAC * 100) + 'pct' +
    '_maxCloudFP' + Math.round(MAX_CLOUD_FRAC_IN_FOOTPRINT * 100) + 'pct' +
    '_R' + WORK_RADIUS_KM + 'km';

Map.addLayer(multiYearMedian, {min: 15, max: 45, palette: ['blue', 'yellow', 'red']}, 'LST_Harmonized_MultiYear', true);

Export.image.toDrive({
  image: multiYearMedian,
  description: multiYearPrefix,
  folder: DRIVE_FOLDER_RASTER,
  fileNamePrefix: multiYearPrefix,
  region: AOI,
  crs: OUT_CRS,
  scale: OUT_SCALE_M,
  maxPixels: MAX_PIXELS
});

