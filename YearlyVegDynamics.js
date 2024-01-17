// import variables including assets/NEONwatersheds as table
// study area is a polygon covering the entire contiguous US
var LandCoverDynamics = ee.ImageCollection("MODIS/061/MCD12Q2");
  
var watersheds = table.select('SiteID');

//calculate the mean value of minimum EVI over 5 years
var EVI_image_min = ee.Image(LandCoverDynamics
        .filterDate('2015-01-01', '2019-12-31')
        .mean().clip(watersheds).select('EVI_Minimum_1')
      );

// Reduce the image to get the mean value for each watershed
var EVImin = EVI_image_min.reduceRegions({
  collection: watersheds,
  reducer: ee.Reducer.mean(),
  scale: 30 // Adjust the scale according to your image resolution
});

// Print the mean values
print(EVImin, 'EVImin');

Export.table.toDrive({
  collection: EVImin,
  description:'min_EVI',
  fileFormat: 'csv',
  selectors: ['SiteID', 'mean']
});


// EVI amplitude
var EVI_image_amp = ee.Image(LandCoverDynamics
        .filterDate('2015-01-01', '2019-12-31')
        .mean().clip(watersheds).select('EVI_Amplitude_1')
      );

var EVIamp = EVI_image_amp.reduceRegions({
  collection: watersheds,
  reducer: ee.Reducer.mean(),
  scale: 30 
});

print(EVIamp, 'EVIamp');

Export.table.toDrive({
  collection: EVIamp,
  description:'amp_EVI',
  fileFormat: 'csv',
  selectors: ['SiteID', 'mean']
});


// EVI area under the curve
var EVI_image_area = ee.Image(LandCoverDynamics
        .filterDate('2015-01-01', '2019-12-31')
        .mean().clip(watersheds).select('EVI_Area_1')
      );

var EVIarea = EVI_image_area.reduceRegions({
  collection: watersheds,
  reducer: ee.Reducer.mean(),
  scale: 30 
});

print(EVIarea, 'EVIarea');

Export.table.toDrive({
  collection: EVIarea,
  description:'area_EVI',
  fileFormat: 'csv',
  selectors: ['SiteID', 'mean']
});


// EVI greenup
var EVI_image_greenup = ee.Image(LandCoverDynamics
        .filterDate('2015-01-01', '2019-12-31')
        .mean().clip(watersheds).select('Greenup_1')
      );

var EVIgreenup = EVI_image_greenup.reduceRegions({
  collection: watersheds,
  reducer: ee.Reducer.mean(),
  scale: 30 
});

print(EVIgreenup, 'EVIgreenup');

Export.table.toDrive({
  collection: EVIgreenup,
  description:'greenup_EVI',
  fileFormat: 'csv',
  selectors: ['SiteID', 'mean']
});

// EVI dormany
var EVI_image_dormancy = ee.Image(LandCoverDynamics
        .filterDate('2015-01-01', '2019-12-31')
        .mean().clip(watersheds).select('Dormancy_1')
      );

var EVIdormancy = EVI_image_dormancy.reduceRegions({
  collection: watersheds,
  reducer: ee.Reducer.mean(),
  scale: 30 
});

print(EVIdormancy, 'EVIdormancy');

Export.table.toDrive({
  collection: EVIdormancy,
  description:'dormancy_EVI',
  fileFormat: 'csv',
  selectors: ['SiteID', 'mean']
});