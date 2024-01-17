var qualityMask = function(im) {
  return im.updateMask(im.select('quality_flag').eq(1))
      .updateMask(im.select('degrade_flag').eq(0));
};
var dataset = ee.ImageCollection('LARSE/GEDI/GEDI02_A_002_MONTHLY')
                  .map(qualityMask)
                  .select('rh98').mosaic();

print('GEDI:', dataset)

var gediVis = {
  min: 1,
  max: 60,
  palette: 'darkred,red,orange,green,darkgreen',
};
Map.setCenter(-74.803466, -9.342209, 10);
Map.addLayer(dataset, gediVis, 'rh98');

var means = dataset.reduceRegions({
    collection: table,
    reducer: ee.Reducer.mean(),
    scale: 30,
  });
  
var stdev = dataset.reduceRegions({
    collection: table,
    reducer: ee.Reducer.stdDev(),
    scale: 30,
  });
  
print('means:', means)
print('standard deviation:', stdev)


Export.table.toDrive({
  collection: means,
  description:'means_GEDI',
  fileFormat: 'csv',
  selectors: ['SiteID', 'mean']
});


Export.table.toDrive({
  collection: stdev,
  description:'stdev_GEDI',
  fileFormat: 'csv',
  selectors: ['SiteID', 'stdDev']
});