// Classification utility functions

var miscUtils = require('users/igor_cnpy/corescam:miscUtils');
var regressionUtils = require('users/igor_cnpy/corescam:regressionUtils');


exports.getCoefsSample = function(regions, params) {
  /* Get regression coefficients for sample locations */
  
  // Get list of years for training period
  var years = ee.List.sequence(params.get('trainDataStart'), params.get('trainDataEnd')).map(function(int) {
    return ee.Date(ee.String(int).slice(0,4))
  })
  // Make a list containing the years, sample region, and param file
  var yearsRegion = years.map(function(year) {
    return ee.List([year, params])
  })
            
  // Make a list of composites  
  var fullList = yearsRegion.map(miscUtils.getMonthlyCompositesTraining)
          
  // Need to turn the list into a composite. For some reason I
  // can only get this to work with an iteration that needs a starting
  // condition.
  var startingCollection = miscUtils.getMonthlyCompositesTraining(
                            ee.List(
                             [ee.Date(ee.String(params.get('trainDataStart'))
                                .slice(0,4))
                                .advance(-1,'year'), 
                              params]))
          
  // Do the iteration to make the collection      
  var inputsFull = ee.ImageCollection(
                    fullList.iterate(
                      miscUtils.concatenateCollections, 
                      startingCollection))

  // Filter to make sure there's actual data
  var inputsTraining = ee.ImageCollection(inputsFull)
                        .filterMetadata('bands','equals',5)
            
  // Get regression coefficients  
  var minRMSE = params.get('minRMSE')
  var regressionNDFIW = regressionUtils.doRegressionTraining(inputsTraining, 'NDFIW', 1, 0, minRMSE)
  var regressionGV = regressionUtils.doRegressionTraining(inputsTraining, 'GVshade', 1, 0, minRMSE)
  var regressionWater = regressionUtils.doRegressionTraining(inputsTraining, 'Water', 1, 0, minRMSE)
  var regressionSoil = regressionUtils.doRegressionTraining(inputsTraining, 'Soil', 1, 0, minRMSE)
  var regressionShade = regressionUtils.doRegressionTraining(inputsTraining, 'Shade', 1, 0, minRMSE)
          
  var coefs = ee.Image(regressionNDFIW).addBands(
                [regressionGV, regressionWater,
                regressionSoil, regressionShade])
  // Extract the coefficient images at sample locations        
  var reducedCoefs = coefs.sampleRegions({
    collection: regions,
    properties: ['label', 'binaryLabel'], // I REMOVED THIS SEPTEMBER 10
    scale: 30,
    tileScale: 4
  })
  return [reducedCoefs]      

}


exports.getTrainedClassifiers = function(sampleName, params) {
  
  // Get trained classifiers for regression information during training
  // period and for for training samples
  
  // To get a binary forest/non-forest probability we need to add a 'binary label'

  // property to the training data. 
  var forestSamples = ee.FeatureCollection(
                        sampleName.filter(
                          ee.Filter.eq(params.get('propertyName'), ee.Number(params.get('forestLabel')))))
                          .map(function(feat){
                            return feat.set('binaryLabel', 1)
                          })
 
  var nfSamples = ee.FeatureCollection(
                        sampleName.filter(
                          ee.Filter.neq(params.get('propertyName'), ee.Number(params.get('forestLabel')))))
                        .map(function(feat){
                            return feat.set('binaryLabel', 0)
                          })
  var binarySamples = ee.FeatureCollection(forestSamples).merge(ee.FeatureCollection(nfSamples))
  
  // Get coefficients and median endmember fractions for samples
  var trainingMapped = ee.List(exports.getCoefsSample(binarySamples, params))
  var trainingWithCoefs = trainingMapped.get(0)

  // Band names with coefficients
  var trainingBandsCoefs =  ee.Feature(
                             ee.FeatureCollection(trainingWithCoefs)
                             .toList(2)
                             .get(0))
                             .propertyNames()
                             .remove(params.get('propertyName'))
                             .remove('system:index')
                             .remove('binaryLabel')

  // Band names with median
  var trainingBandsMedian = ['GVshade','NDFIW','Water','Shade','Soil']

  // Get trained classifier with the coefs and median data
  var trainedClassifierCoefs = ee.Classifier.smileRandomForest({
                                 numberOfTrees:200})
                                 .train(trainingWithCoefs, params.get('propertyName'), trainingBandsCoefs)
  
  // Binary is for forest/non-forest probability                              
  var trainedClassifierBinary = ee.Classifier.smileRandomForest({
                                  numberOfTrees:200}).setOutputMode('PROBABILITY')
                                  .train(trainingWithCoefs, 'binaryLabel', trainingBandsCoefs)


  return [trainedClassifierCoefs, trainedClassifierBinary]
  
}

