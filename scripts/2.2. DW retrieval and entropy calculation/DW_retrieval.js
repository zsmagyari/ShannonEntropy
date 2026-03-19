var CITY_KEY = 'BUCURESTI'; 

var CITIES = {
  'BUCURESTI':   {name:'București',   lon:26.1025, lat:44.4268, zoom:10},
  'CLUJ_NAPOCA': {name:'Cluj-Napoca', lon:23.5899, lat:46.7712, zoom:11},
  'IASI':        {name:'Iași',        lon:27.5879, lat:47.1585, zoom:11},
  'CONSTANTA':   {name:'Constanța',   lon:28.6348, lat:44.1598, zoom:11},
  'TIMISOARA':   {name:'Timișoara',   lon:21.2287, lat:45.7489, zoom:11},
  'BRASOV':      {name:'Brașov',      lon:25.6012, lat:45.6579, zoom:11},
  'CRAIOVA':     {name:'Craiova',     lon:23.7949, lat:44.3302, zoom:11}
};

var city = CITIES[CITY_KEY];
if (!city) throw new Error('Unknown CITY_KEY: ' + CITY_KEY);

var WORK_RADIUS_KM = 35; 
var AOI = ee.Geometry.Point([city.lon, city.lat]).buffer(WORK_RADIUS_KM * 1000);

Map.setCenter(city.lon, city.lat, city.zoom);
Map.addLayer(AOI, {}, 'AOI');

var YEARS = [2021, 2022, 2023, 2024, 2025];
var MONTH_START = 6;   // JJA
var MONTH_END   = 8;

var PROB_BANDS = [
  'water','trees','grass','flooded_vegetation','crops',
  'shrub_and_scrub','built','bare','snow_and_ice'
];

var COMPOSITE = 'median';

var DRIVE_FOLDER = 'GEE_DW_yearly_JJA_DWmaskOnly';
var OUT_CRS      = 'EPSG:3035';
var OUT_SCALE_M  = 10;
var MAX_PIXELS   = 1e13;
var FILL_NODATA  = false;
var NODATA_VALUE = -9999;

var DW = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1')
  .filterBounds(AOI)
  .select(PROB_BANDS);

var annualValidImages = [];

YEARS.forEach(function(y) {
  var start = ee.Date.fromYMD(y, 1, 1);
  var end   = ee.Date.fromYMD(y + 1, 1, 1);

  var dwY = DW
    .filterDate(start, end)
    .filter(ee.Filter.calendarRange(MONTH_START, MONTH_END, 'month'));

  print('Year', y, 'DW scenes (JJA):', dwY.size());
  
  var dwKeyed = dwY.map(function(img){
    return ee.Image(img).set('dateKey', ee.Date(img.get('system:time_start')).format('YYYYMMdd'));
  });

  var keys = ee.List(dwKeyed.aggregate_array('dateKey')).distinct().sort();

  var dailyValid = ee.ImageCollection.fromImages(
    keys.map(function(k){
      k = ee.String(k);
      var day = dwKeyed.filter(ee.Filter.eq('dateKey', k));
      var v = day.map(function(img){
          return img.select('built').mask().rename('valid').unmask(0);
        }).max().rename('valid');
      return v.set('dateKey', k);
    })
  );

  var nValid = dailyValid.sum().rename('N_valid').clip(AOI);
  
  annualValidImages.push(nValid);

  var annual;
  if (COMPOSITE === 'mean') {
    annual = dwY.mean();
  } else {
    annual = dwY.median();
  }

  annual = annual.select(PROB_BANDS).clip(AOI);
  if (FILL_NODATA) {
    annual = annual.unmask(NODATA_VALUE);
  }

  Map.addLayer(annual.select('built'), {min:0, max:1}, 'built_' + y, false);
});

var multiYearNValid = ee.ImageCollection(annualValidImages).sum().rename('N_valid');
var startMulti = ee.Date.fromYMD(YEARS[0], 1, 1);
var endMulti   = ee.Date.fromYMD(YEARS[YEARS.length-1] + 1, 1, 1);

var dwMulti = DW
  .filterDate(startMulti, endMulti)
  .filter(ee.Filter.calendarRange(MONTH_START, MONTH_END, 'month'));

print('Multi-Year Total Scenes:', dwMulti.size());

var multiAnnual;
if (COMPOSITE === 'mean') {
  multiAnnual = dwMulti.mean();
} else {
  multiAnnual = dwMulti.median();
}

multiAnnual = multiAnnual.select(PROB_BANDS).clip(AOI);

if (FILL_NODATA) {
  multiAnnual = multiAnnual.unmask(NODATA_VALUE);
}

Map.addLayer(multiAnnual.select('built'), {min:0, max:1}, 'MULTIYEAR_built', true);

var multiPrefix = 'RO_' + CITY_KEY + '_DW_JJA_MULTIYEAR_' + YEARS[0] + '-' + YEARS[YEARS.length-1] +
                  '_' + COMPOSITE +
                  '_R' + WORK_RADIUS_KM + 'km' +
                  (FILL_NODATA ? '_NODATA' : '_MASK');

Export.image.toDrive({
  image: multiYearNValid,
  description: multiPrefix + '_N_valid',
  folder: DRIVE_FOLDER,
  fileNamePrefix: multiPrefix + '_N_valid',
  region: AOI,
  crs: OUT_CRS,
  scale: OUT_SCALE_M,
  maxPixels: MAX_PIXELS
});

PROB_BANDS.forEach(function(b) {
  Export.image.toDrive({
    image: multiAnnual.select(b),
    description: multiPrefix + '_' + b,
    folder: DRIVE_FOLDER,
    fileNamePrefix: multiPrefix + '_' + b,
    region: AOI,
    crs: OUT_CRS,
    scale: OUT_SCALE_M,
    maxPixels: MAX_PIXELS
  });
});
