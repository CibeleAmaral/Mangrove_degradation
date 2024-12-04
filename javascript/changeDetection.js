// Continuous Degradation Detection (CODED)
// Version 0.2

var miscUtils = require('users/igor_cnpy/corescam:miscUtils');
var regressionUtils = require('users/igor_cnpy/corescam:regressionUtils');
var changeUtils = require('users/igor_cnpy/corescam:changeUtils');
var class_utils = require('users/igor_cnpy/corescam:classUtils');


var coded = function(inputsTraining, inputsMonitoring, pgeo, 
                     params, coefClassifier, binaryClassifier) {
  
  /*  Run regression and change detection for CODED  */
  
  // Training and calculating residuals
  
  // regressionData = regression coefficients and observation residuals
  var regressionData = ee.List(regressionUtils.allRegression(inputsTraining, 
                                                             inputsMonitoring, 
                                                             1, 
                                                             0, 
                                                             .015))
  
  // residaulsNDFIW = NDFIW residuals during monitoring period based on 
  // regression from training period
  var residualsNDFIW = ee.ImageCollection(regressionData.get(0))
  
  // allCoefs = NDFIW, shade, soil, GV, and Water coefficients for classification
  var allCoefs = ee.Image(regressionData.get(1))
  
  // rmses = NDFIW, shade, soil, GV, and Water RMSE for classification
  var rmses = ee.Image(regressionData.get(2))
  
  // residualsNDFIWCu = NDFIW residuals based on regression fit from monitoring period,
  // used for CUSUM change detection
  var residualsNDFIWCu = ee.ImageCollection(regressionData.get(3))


  // Classifification
  
  // classifiedTraining = land cover classification during training period
  var classifiedTraining = allCoefs.classify(coefClassifier)
  
  // forest = forest flag (1 = forest)
  var forest = classifiedTraining.eq(ee.Image.constant(params.get('forestLabel'))) 
  
  // forestProbabilityTrain = binary forest classification probability
  var forestProbabilityTrain = allCoefs.classify(binaryClassifier)


  // Change detection
  
  // collectionIterate = input data combined with residuals
  var collectionIterate = inputsMonitoring.map(function(image) {
    var theDate = image.date()
    var dateBefore = theDate.advance(-1, 'day')
    var dateAfter = theDate.advance(1,'day')
    var regressionImage = residualsNDFIW.filterDate(dateBefore, dateAfter).first()
    var regressionImageCu = residualsNDFIWCu.filterDate(dateBefore, dateAfter).first()
    return image.addBands([ee.Image(regressionImage).updateMask(forest), 
                           ee.Image(regressionImageCu).rename('residualCUSUM')])
  })
  
  // statusList = empty image to hold change information during monitor iteration
  var statusList = [ee.Image(0).addBands([ee.Image(0),ee.Image(0),ee.Image(0)])
                      .rename([
                        'consecutive',
                        'changeFlag',
                        'magnitude',
                        'changeDate'])
                      .cast({
                        'consecutive': 'short',
                        'changeFlag': 'short',
                        'magnitude': 'float',
                        'changeDate': 'float'
                      })
                      .setMulti({
                        'system:time_start': 0,
                        'consec': ee.Number(params.get('consec')).short(),
                        'thresh': ee.Number(params.get('thresh')).short()})]

  // results = change detection results     
  var results = ee.ImageCollection(
                    ee.List(collectionIterate.select(['normalizedRes'])
                    .iterate(changeUtils.monitorFunction,statusList)))
  return [results, classifiedTraining, allCoefs, forestProbabilityTrain]

  //return [results, forest, allCoefs, forestProbabilityTrain]
}



var mappingFunction = function(years, mappingYears, iCol, length, region, params,
                               coefClassifier, binaryClassifier, skipLength) {
  /* Function to handle mapping over every year in study period */
  // yearlyData = algorithm outputs for every year in study period
  var yearlyData = mappingYears.map(function(year) {
    
    // Inputs
    
    // trainEnd = last year in training period
    var trainEnd = ee.Number(year).add(ee.Number(length)).toInt()
    
    // monitorEnd = last year in monitor period, 2 years after trainEnd
    var monitorEnd = trainEnd.add(ee.Number(params.get('window'))).toInt() //TODO
    
    // trainBegDate = beginning date of training period
    var trainBegDate = ee.Date(ee.String(ee.Number(year).toInt()).cat('-01-01'))
    
    // trainEndDate = end date of training period
    var trainEndDate = ee.Date(ee.String(trainEnd).cat('-01-01'))
    // monitorEndDate = end date of monitoring period
    var monitorEndDate = ee.Date(ee.String(monitorEnd).cat('-12-31'))
    // trainingDataFull = inputs filtered by training period
    var trainingDataFull = iCol.filterDate(trainBegDate, trainEndDate)
    
    // monitorDataFull = inputs filtered by monitoring period
    var monitorDataFull = iCol.filterDate(trainEndDate, monitorEndDate)

    // maskedObservations = inputs masked if there is not enough training data
    var maskedObservations = ee.List(miscUtils.maskLowObs(trainingDataFull, 
                                                          monitorDataFull,
                                                          params))
    
    // trainingData = training data masked if there was not enough observations
    var trainingData = ee.ImageCollection(maskedObservations.get(0))
    
    // monitorData = monitoring data masked if there was not enough observations
    var monitorData = ee.ImageCollection(maskedObservations.get(1))


    // Change Detection

    // results = change detection results
    var results = ee.List(coded(trainingData, monitorData, 
                              region, params, coefClassifier,
                              binaryClassifier))
                              
    // resultList = change detection results as list
    var resultList = ee.List(miscUtils.cleanupResults(results))

    // resultImage = change detection results as image
    var resultImage = ee.Image(resultList.get(0)).unmask()
    
    // changeDate = date of change image
    var changeDate = ee.Image(resultList.get(1)).unmask()
    
    // forestMask = forest flag (1=forest) during training period
    var forestMask = ee.Image(resultList.get(2)).unmask()
    
    
    // Retrain
    
    // inputsRetrain = inputs to use for retrain period
    var inputsRetrain = ee.List(miscUtils.GetInputsRetrain(iCol, 
                          changeDate, 
                          2, 
                          params, length, skipLength))
  
    // trainingRetrain = training inputs post disturbance
    var trainingRetrain = inputsRetrain.get(0)
    
    // monitorRetrain = monitor inputs postchange
    var monitoringRetrain = inputsRetrain.get(1)
  
    // regressionRetrain = Regression coefficients for retrian period
    var regressionRetrain = ee.List(
                              regressionUtils.allRegression(
                                trainingRetrain, 
                                monitoringRetrain, 
                                1, // number harmonics
                                0, // no trend
                                .015))
                                
    // retrainCoefs = regression coefficients post disturbance
    var retrainCoefs = ee.Image(regressionRetrain.get(1))
    
    // classifiedRetrain = post disturbance classification
    var classifiedRetrain = retrainCoefs.classify(coefClassifier)
                              .rename('classAfter')
    
    // trainNDFIW = NDFIW during pre-disturbance training period
    var trainNDFIW = ee.Image(results.get(2))
    
    // difNDFIW = difference in NDFIW before and after disturbance
    var difNDFIW = retrainCoefs.select('coef_constant_NDFIW')
                    .divide(trainNDFIW.select('coef_constant_NDFIW'))
                    .multiply(ee.Image(100))
                    .rename('difNDFIW')
    
    // return change results, post-disturbance classification, difference in 
    // NDFIW, and the forest mask
    return resultImage.addBands([classifiedRetrain, difNDFIW, forestMask])
  })

  var outputImages = ee.List(miscUtils.organizeOutputImages(yearlyData,
                                                            mappingYears,
                                                            years,
                                                            params))
  
  // Finally convert events labeled as deg but then not classified as forest 
  // for many years
  var postDistFinal = miscUtils.convertFalseDeg(ee.Image(outputImages.get(2)),
                                                ee.Image(outputImages.get(4)), 
                                                params)

                                                            
  // Convert output to array image
    // var changeOnly1 = miscUtils.collapseList(yearlyData, 'changeDate')
  // var changeOnly2 = ee.Image(outputImages)

  // var arr1 = ee.Image(changeOnly).unmask().toArray().toArray(1)
  // return [changeOnly1, changeOnly2]
  
  // arr1 = change flag as array image

  var arr1 = ee.Image(outputImages.get(0)).unmask().toArray().toArray(1)

  
  // arr2 = change magnitude as array image
  var arr2 = ee.Image(outputImages.get(1)).unmask().toArray().toArray(1)
  
  // arr3 = post disturbance class as array image
  var arr3 = ee.Image(postDistFinal).unmask().toArray().toArray(1)
  
  // arr4 = difference in NDFIW as array image
  var arr4 = ee.Image(outputImages.get(3)).unmask().toArray().toArray(1)
  
  // arr5 = forest flag as array image
  var arr5 = ee.Image(outputImages.get(4)).unmask().toArray().toArray(1)
  
  var arrayOutput = arr1.arrayCat(arr2, 1)
                        .arrayCat(arr3, 1)
                        .arrayCat(arr4, 1)
                        .arrayCat(arr5, 1)
                        .toShort()

  return arrayOutput
  
}  



exports.submitCODED = function(region, params, samples) {

  /* Function to manage the running of CDD across iterations through the data */
  

  // Parameters
  
  // start = first year of training
  var start = ee.Number(params.get('start'))
  
  // end = last year of monitoring
  var end = ee.Number(params.get('end'))
  
  // trainLength = number of years for training
  var trainLength = ee.Number(params.get('trainLength'))
  
  // lastYear = last year to begin training
  var lastYear = ee.Algorithms.If((end.subtract(trainLength).subtract(2)).gt(0),
                                  // end.subtract(trainLength).subtract(1), //1
                                  end.subtract(trainLength), //1
                                   end) //0
                                   
  var skipLength = ee.Algorithms.If(params.contains('skipLength'),
    ee.Number(params.get('skipLength')),
    2)
    skipLength = ee.Number(skipLength)
  // Inputs
  
  // Training parameter name
  params = ee.Algorithms.If(params.contains('propertyName'),
    params,
    params.set('propertyName','label'))

  params = ee.Dictionary(params)
  // years = full list of years including the end of the time series
  var years = ee.List.sequence(start, end)

  // mappingYears = list year of years that can begin with training
  var mappingYears = ee.List.sequence(start, lastYear, 2)

  // yearsRegion = list of years with study region and params attached
  var yearsRegion = ee.List(miscUtils.getYears(start, end.add(1), region, params))

  // fullList = composites of every month in study period
  var fullList = yearsRegion.map(miscUtils.getMonthlyComposites)

  // startingCollection = temporary collection to merge fullList, which is a list
  var startingCollection = miscUtils.getMonthlyComposites(ee.List([ee.Date('1970'), region, params]))

  // inputsFull = fullList converted to an image collection using an iteration
  var inputsFull = ee.ImageCollection(
                    fullList.iterate(
                      miscUtils.concatenateCollections, 
                      startingCollection))
  // inputs = only keep monthly composites where there is data
  var inputs = ee.ImageCollection(inputsFull)
          .filterMetadata('bands','equals',5)

  // return mappingYears}
  // bothClassifiers = landcover classifiers
  var bothClassifier = ee.List(class_utils.getTrainedClassifiers(samples, params))
  // coefClassifier = primary classifier used to classify coefficients into land covers
  var coefClassifier = bothClassifier.get(0)

  // binaryClassifier = 0-1 classifier to get forest probability
  var binaryClassifier = bothClassifier.get(1)

  // results = output of algorithm from mapping over every year to look for change
  var results = mappingFunction(years, mappingYears, inputs, trainLength, region, params,
                                coefClassifier, binaryClassifier, skipLength)
  return results
}
