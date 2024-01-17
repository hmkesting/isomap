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

var histograms = dataset.reduceRegions({
    collection: table,
    reducer: ee.Reducer.fixedHistogram(0, 115, 23),
    scale: 30,
  });
  
print('histograms:', histograms)

var histogram_values = histograms.map(function(feat){
  feat = ee.Feature(feat);
  // cast to an array and get the bucket means and counts
  var hist = ee.Array(feat.get('histogram'));
  var means = hist.slice(1, 1, 2).project([0]);
  var counts = hist.slice(1, 0, 1).project([0]); // eventually cast to a list using toList()
  return feat.set('means', means, 'counts', counts);
});

Export.table.toDrive({
  collection: histogram_values,
  description:'hist_GEDI',
  fileFormat: 'csv',
  selectors: ['SiteID', 'counts', 'histogram', 'means']
});