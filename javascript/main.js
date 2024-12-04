// It is necessary to declare new imports as points, with the classes of: mangroves, water, soil and salt marsh

// Region of Interest 
var roi = ee.FeatureCollection('path/to/asset');

// Parameters
var params = ee.Dictionary({
    'cfThreshold': .05,
    'consec': 3,
    'thresh': 3,
    'start': 2000,
    'end': 2022,
    'trainDataEnd': 1999,
    'trainDataStart': 1997,
    'trainLength': 3,
    'forestLabel': 1,
    'window': 2,
    'minYears': 3,
    'minRMSE': .015,
    'numChanges': 6,
    'minObs': 5
})

// required scripts sheets
var codedUtils = require('users/igor_cnpy/corescam:changeDetection')
var dataUtils = require('users/igor_cnpy/corescam:dataUtils')


// Training data
var trainingData = mangrove.merge(soil).merge(water).merge(saltMarsh)


// Output names
var assetName = 'Mangrove_Degradation_for_'+params.get('start').getInfo()+'_'+params.get('end').getInfo()

// Monitoring main function
var results = codedUtils.submitCODED(roi, params, trainingData)


// Turn array columns into images
var disturbances = dataUtils.makeImage(results, 0, 'dist_', params.get('start'), params.get('end'))
var magnitude = dataUtils.makeImage(results, 1, 'mag_', params.get('start'), params.get('end'))
var postChange = dataUtils.makeImage(results, 2, 'post_', params.get('start'), params.get('end'))
var difference = dataUtils.makeImage(results, 3, 'dif_', params.get('start'), params.get('end'))
var forestFlag = dataUtils.makeImage(results, 4, 'forest_', params.get('start'), params.get('end'))

var disturbanceBands = disturbances.addBands([magnitude, postChange, difference])

var save_output = ee.Image(dataUtils.reduceBands(ee.Image(disturbanceBands), params)
                           .addBands(forestFlag.select(0)) // Forest flag for first year
                           .setMulti(params))



// Save the results as assets
Export.image.toAsset({
  image: save_output,
  description: assetName,
  assetId: assetName,
  maxPixels: 1e13,
  scale: 30,
  region: roi,
  pyramidingPolicy: {
  '.default': 'mode'
  }
})


 
