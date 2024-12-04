// Miscellaneous utily functions for CODED

exports.convertFalseDeg = function(postDistBands, forestBands, params) {
  /* Convert locations identified as degradation but not classified as 
  forest for many years after the disturbance */
  
  // Where it starts as a forest
  var startForest = forestBands.select(0)
                      .eq(ee.Image(ee.Number(params.get('forestLabel'))))
  
  // Where there is only one change
  var singleChange = postDistBands.gt(0).reduce(ee.Reducer.sum()).eq(1)
  
  // Where there was deg and only one change
  var oneChangeDeg = singleChange.eq(1)
                      .and(postDistBands.reduce(ee.Reducer.max())
                        .eq(ee.Image(ee.Number(params.get('forestLabel')))))
  
  // Where it was classified as nonforest for trainLength + 2
//  var multipleNF = forestBands.eq(0)
  var multipleNF = forestBands.neq(ee.Image(ee.Number(params.get('forestLabel'))))
                    .reduce(ee.Reducer.sum())
                    .gte(ee.Image(ee.Number(params.get('trainLength')).add(2)))
                    .and(oneChangeDeg.eq(1))
  
                    
  var fixedClasses = postDistBands
                        .multiply(multipleNF)
                        .multiply(ee.Image(ee.Number(params.get('forestLabel')))) 
                        .add(
                          postDistBands
                          .multiply(multipleNF.eq(0)))
                          
  return fixedClasses
}

exports.organizeOutputImages = function(yearlyData,mappingYears, years, params) {
    // Results from each year need to be dissolved into something logical
  var changeDates = exports.collapseList(yearlyData, 'changeDate')
  var magnitude = exports.collapseList(yearlyData, 'magnitude')
  var difNDFIW = exports.collapseList(yearlyData, 'difNDFIW')
  var classAfter = exports.collapseList(yearlyData, 'classAfter')
  var forest = exports.collapseList(yearlyData, 'forestTraining')
  
  // Get forest value for every year
  var forestNeeded = years.length().subtract(mappingYears.length()).subtract(1)
  var newBands = ee.List.sequence(0, forestNeeded).map(
                  function(i) {
                    return ee.Image(99).rename('extra')
                  })
                  
  var newBandsImage = exports.collapseList(newBands, 'extra')                
  var forestFull = forest.addBands(newBandsImage)
  
  // remove changes higher in frequency that minYears

  var bandList = ee.List.sequence(0, years.length())
  
  // fix the bands to align with change dates // TODO
  var dataTimeFixed = years.map(function(i) {
    // if there is a change this year
       
    var changeThisYear = changeDates.eq(ee.Number(i).toInt())
                              .reduce(
                                ee.Reducer.max())
                              .rename('disturbance')
    
    var changeMask = changeDates.eq(ee.Number(i).toInt())
                              
    var magnitudeThisYear = ee.Image(magnitude.multiply(changeMask))
                              .reduce(
                                ee.Reducer.max())
                              .rename('magnitude')
                              
    var classThisYear = classAfter.multiply(changeMask)
                              .reduce(
                                ee.Reducer.max())
                              .rename('class')
                              
    var difThisYear = difNDFIW.multiply(changeMask)
                              .reduce(
                                ee.Reducer.max())
                              .rename('dif')

        
    return changeThisYear.addBands([magnitudeThisYear,classThisYear, difThisYear])

  })

  // Get years between disturbances
  var minYears = ee.Number(params.get('minYears'))
  var startImage = ee.List([ee.Image(0)])
  var iterList = exports.cumulativeDists(startImage, dataTimeFixed, minYears, 'disturbance')
  var iterImage = exports.collapseList(iterList, 'changeDate')
  
  // Now turn to zero anywhere that's beyond minYears
  var numYears = ee.List.sequence(0, iterList.length().subtract(1)) //TODO
  // var numYears = ee.List.sequence(0, iterList.length())

  var changes = exports.removeExtras(iterList, dataTimeFixed, numYears, minYears, 'disturbance')
  var mag = exports.removeExtras(iterList, dataTimeFixed, numYears, minYears, 'magnitude')
  var classes = exports.removeExtras(iterList, dataTimeFixed, numYears, minYears, 'class')
  var dif = exports.removeExtras(iterList, dataTimeFixed, numYears, minYears, 'dif')
  
  var changeAllFixed = exports.collapseList(changes, 0)
  var magAllFixed = exports.collapseList(mag, 0)
  var classAllFixed = exports.collapseList(classes, 0)
  var difAllFixed = exports.collapseList(dif, 0)
  
  var lowDifChange = difAllFixed.gte(80)
  var highDifChange = difAllFixed.lt(80)
                        .and(changeAllFixed.eq(1))
                        
  var deg = lowDifChange.eq(1)
  var def = highDifChange.eq(1)
  var classConverted = deg.multiply(ee.Image(ee.Number(params.get('forestLabel'))))
                        .add(classAllFixed.multiply(def))
                        
  // Finally, remove changes if only F for one year
  // Where it starts as non-forest
  var startNonForest = forestFull.select(0)
                      .neq(ee.Image(ee.Number(params.get('forestLabel'))))
  
  // Where there is only one change
  //var singleChange = classConverted.gt(0).reduce(ee.Reducer.sum()).eq(1)
  var singleChange = changeAllFixed.gt(0).reduce(ee.Reducer.sum()).eq(1)

  
  // Where it was classified as forest for only one period
  var singleF = forestFull.eq(ee.Image(ee.Number(params.get('forestLabel'))))
                    .reduce(ee.Reducer.sum())
                    .eq(ee.Image(1))
                    .and(singleChange.eq(1))
                    .and(startNonForest.eq(1))

  var classFinal = classConverted
                        .multiply(singleF)
                        .multiply(ee.Image(0))
                        .add(classConverted.multiply(singleF.eq(0)))
                        
  var changeFinal = changeAllFixed
                        .multiply(singleF)
                        .multiply(ee.Image(0))
                        .add(changeAllFixed.multiply(singleF.eq(0)))
                        
  var magFinal = magAllFixed
                        .multiply(singleF)
                        .multiply(ee.Image(0))
                        .add(magAllFixed.multiply(singleF.eq(0)))
                        
  var difFinal = difAllFixed
                        .multiply(singleF)
                        .multiply(ee.Image(0))
                        .add(difAllFixed.multiply(singleF.eq(0)))

  return [changeFinal, magFinal, classFinal, difFinal, forestFull]
}
                                                          

exports.cleanupResults = function(list) {
  /* Turn the results into something that makes more sense */

  var results = ee.List(list)
  var resultsColTemp = ee.ImageCollection(results.get(0))
  var forestMask = ee.Image(results.get(1)).rename('forestTraining')
  
  var resultsMax = resultsColTemp.max()

  // Where change was detected
  var changeFlag = resultsMax.select('changeFlag')
  
  // Hack to get the date the change was first detected, it is only recorded
  // if there was no previous change and on first observation passed threshold
  var changeDate = resultsMax.select('changeDate')
                    .add(1970).toInt() 
                    .multiply(changeFlag)

  // Add true change date band instead of years since 1970
  var resultsCol = resultsColTemp.map(function(image) {
                     return image.addBands(changeDate.rename('trueChangeDate'))
                   })

  // Change magnitude
  var changeMagnitude = ee.Image(resultsCol.map(exports.maskMagnitude)
                          .mean())
                          .select('magnitude')
                          .multiply(changeFlag) 
                          
  var resultImage = changeDate.addBands(changeMagnitude)
                      .rename([
                        'changeDate',
                        'magnitude'])
                        
  return [resultImage, changeDate, forestMask]
}


exports.maskLowObs = function(trainCol, monCol, params) {
  /* Mask observations with too few training data */
  var imageCount = ee.ImageCollection(trainCol).count()
  var imageCount2 = ee.ImageCollection(monCol).count()
  // var obsMask = imageCount.select('NDFIW').gt(12)
  // minObs
  var obsMask = imageCount.select('NDFIW')
                  .gte(ee.Image(ee.Number(params.get('minObs'))))
  var obsMaskMonitor = imageCount2.select('NDFIW')
                  .gte(ee.Image(ee.Number(params.get('minObs'))))
                  
  var maskedTrain = trainCol.map(function(i) {
                  return ee.Image(i).updateMask(obsMask)
  })
  
  var maskedMon = monCol.map(function(i) {
                  return ee.Image(i).updateMask(obsMaskMonitor)
  })
  return [maskedTrain, maskedMon]
}

exports.cumulativeDists = function(startImage, imList, minYears, band) {
  /* Get years since last disturbance */
  var iterFunc = function(curImage, imageList) {
    // Get last image
    var lastImage = ee.List(imageList).get(-1)
    
    // Current image
    var curChange = ee.Image(curImage).select(band).gt(0)
                      .or(ee.Image(lastImage).gt(0))
    
    // Get cumulative
    var acceptable = ee.Image(lastImage).lt(minYears)
                      
    var cumuChange = ee.Image(lastImage).add(curChange)
                      .multiply(acceptable) // only keep if currently change
                      .rename('changeDate')
                      
    return ee.List(imageList).add(cumuChange)

  }
  var iterListFull = imList.iterate(iterFunc, startImage)
  var iterList = ee.List(iterListFull).remove(ee.List(iterListFull).get(0))
  // var iterList = ee.List(iterListFull)

  return iterList
}

exports.revertList = function(im) {
  /* Convert and image to a list of images with each image being a band
  from the original image */
  var bn = ee.Image(im).bandNames()
  var imageList = bn.map(function(i) {
    return ee.Image(ee.Image(im).select([i]))
  })
  return ee.List(imageList)
}

exports.removeExtras = function(changeList, toFixList, years, minYears, band) {
  /* Remove changes that fall beyond minYears apart from one another */
  var removed = ee.List(years).map(function(i) {
    var curImage = ee.Image(ee.List(toFixList).get(i)).select(band)
    var changeImage = ee.Image(ee.List(changeList).get(i)) 
    
    var acceptable = changeImage.eq(1) 
    return curImage.multiply(acceptable)
    
  })
  return removed
}

exports.renameBands = function(image, prefix, years){
  /* Rename the bands of an image */
  var newNames = years.map(function(i) {
    return ee.String(prefix).cat(ee.String(ee.Number(i).toInt()))
  })
  return ee.Image(image).rename(newNames)
}

exports.collapseList = function(list, band) {
  /* Convert a list to an image with multiple bands */
  var firstImage = ee.Image(0)
  var collapsed = list.iterate(function(newImage, outImage) {
    return ee.Image(outImage).addBands(ee.Image(newImage).select(band))
    },
    firstImage)
  return ee.Image(collapsed).select(ee.Image(collapsed).bandNames().remove('constant'))
}

exports.GetInputsRetrain = function(iCol, dateImage, iteration, params, trainLength, skipLength) {
  /* Get inputs for classification after a change */
  var hasChanged = ee.Image(dateImage).gt(0)
                      
  var changeDateReal = ee.Image(dateImage)//.add(1970)
  var trainPeriodBeg = changeDateReal.add(ee.Image(ee.Number(skipLength)))
  var trainPeriodEnd = changeDateReal.add(trainLength)
  
  var newTraining = ee.ImageCollection(iCol).map(function(image) {
    var theDate = image.metadata('system:time_start')
                   .divide(ee.Image(315576e5)).add(1970)
    var maskDate = theDate.gt(trainPeriodBeg)
                    .and(theDate.lt(trainPeriodEnd))
                    .and(hasChanged.eq(1))
    return image.updateMask(maskDate)
  })
  
  var newMonitoring = ee.ImageCollection(iCol).map(function(image) {
    var theDate = image.metadata('system:time_start')
                   .divide(ee.Image(315576e5)).add(1970)
    var maskDate = theDate.gt(trainPeriodEnd)
                    .and(hasChanged.eq(1))
    return image.updateMask(maskDate)
  })
  
  return [newTraining, newMonitoring]
}

exports.maskMagnitude = function(image) {
  /* Change magnitude is always recorded, so mask out areas of no change */
  var changeDateName = 'trueChangeDate'
  var changeDate = ee.Image(image).select(changeDateName)
  var hasChanged = changeDate.gt(0)
        
  var imageDate = image.metadata('system:time_start')
                   .divide(ee.Image(315576e5))
                   .add(ee.Image(1970))

  return image.updateMask(imageDate.gte(changeDate)
                          .and(hasChanged.eq(1))
                          .and(image.neq(0)))
                          .abs()
}


exports.makeVariables = function(image) {
  /*  Computes the predictors and the response from the input.  */
  var year = ee.Image(image.date().difference(ee.Date('1970-01-01'), 'year'))
  var season = year.multiply(2 * Math.PI)
  return image.select().addBands(ee.Image(1)).addBands(
    season.sin().rename(['sin'])).addBands(
    season.cos().rename(['cos'])).addBands(
    image).toFloat()
}


exports.getInputsRegion = function(start, end, region, params) {
  /* Get inputs for a region of interest */

  var roi = region;

    var SD = start;
  var ED = end;

  // bandas landsat 5
  var bands_5 = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'];
  // bandas landsat 8
  var bands_8 = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']; 
  
  // Endmembers LT05 e LC08

  var cloud_5 = [16022.125,	7781.473195,	7843.746653,	8623.887614,	6047.753435,	3765.897382]
  var gv_5 = [261.5757123,	577.9617759,	366.6097816,	5291.62735,	2121.13877,	770.4923314]
  var soil_5 = [1036.654025,	1466.611126,	1777.302094,	2795.225033,	4631.590805,	3619.778338]
  var wat_5 = [2090.036621,	2955.367136,	1767.187081,	232.7541696,	35.7418541,	21.8578176]
  var shade_5 = [0, 0, 0, 0, 0, 0];
  
  var cloud_7 = [16022.125,	7781.473195,	7843.746653,	8623.887614,	6047.753435,	3765.897382]
  var gv_7 = [261.5757123,	577.9617759,	366.6097816,	5291.62735,	2121.13877,	770.4923314]
  var soil_7 = [1036.654025,	1466.611126,	1777.302094,	2795.225033,	4631.590805,	3619.778338]
  var wat_7 = [2090.036621,	2955.367136,	1767.187081,	232.7541696,	35.7418541,	21.8578176]
  var shade_7 = [0, 0, 0, 0, 0, 0];
  
  var cloud_8 = [ 6375.76141,	6693.028614,	6780.477957,	6968.257616,	7492.911619,	5969.836111,	4662.671087];
  var gv_8 = [281.2302661,	281.1218856,	645.1652039,	269.2813137,	6280.259655,	2007.282559,	735.9335787];
  var soil_8 = [ 1633.97976,	1923.408204,	2963.690719,	3939.457246,	5013.318084,	6035.979192,	4142.394072];
  var wat_8 = [1118.224128,	1611.110856,	2307.231735,	1311.763986,	113.4761094,	181.8146989,	145.5534865];
  var shade_8 = [0, 0, 0, 0, 0, 0, 0];
  
  
  // filtro de nuvem cloudMaskL457
var cloudMask_L57 = function(image) {
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
  var dilated_cloud = (1 << 1);
  var cloud = (1 << 3);
  var cloud_shadow = (1 << 4);
  var snow = (1 << 5);
  var water = (1 << 7);
  var cloud_conf = (1 << 9);
  var cloud_shadow_conf = (1 << 11);
  var snow_ice_conf = (1 << 12);
  var cirrus_conf = (1 << 15);
  
  // Get the pixel QA band.
  var qa = image.select('QA_PIXEL').rename('clouds');
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloud)
  .and(qa.bitwiseAnd(cloud_conf))
  .or((qa.bitwiseAnd(cloud_shadow))
  .and(qa.bitwiseAnd(cloud_shadow_conf)))
  .or((qa.bitwiseAnd(snow))
  .and(qa.bitwiseAnd(snow_ice_conf)))
  .or(qa.bitwiseAnd(cirrus_conf))
  .or(qa.bitwiseAnd(dilated_cloud))
  .or(qa.bitwiseAnd(water))
  var mask2a = image.mask().reduce(ee.Reducer.min())
    return image.updateMask(mask.not()).updateMask(mask2a);
     
};

// filtro de nuvem cloudMaskL457
var cloudMask_L89 = function(image) {
    // Bits 3 and 5 are cloud shadow and cloud, respectively.
    var dilated_cloud = (1 << 1);
    var cirrus = (1 << 2);
    var cloud = (1 << 3);
    var cloud_shadow = (1 << 4);
    var snow = (1 << 5);
    var water = (1 << 7);
    var cloud_conf = (1 << 9);
    var cloud_shadow_conf = (1 << 11);
    var snow_ice_conf  = (1 << 12);
    var cirrus_conf = (1 << 15);
    
    // Get the pixel QA band.
    var qa = image.select('QA_PIXEL');
    // Both flags should be set to zero, indicating clear conditions.
    var mask = qa.bitwiseAnd(cloud)
                  .and(qa.bitwiseAnd(cloud_conf))
                  .or((qa.bitwiseAnd(cloud_shadow))
                  .and(qa.bitwiseAnd(cloud_shadow_conf)))
                  .or((qa.bitwiseAnd(snow))
                  .and(qa.bitwiseAnd(snow_ice_conf)))
                  .or((qa.bitwiseAnd(cirrus))
                  .and(qa.bitwiseAnd(cirrus_conf)))
                  .or(qa.bitwiseAnd(dilated_cloud))
                  .or(qa.bitwiseAnd(water))
    var mask2a = image.mask().reduce(ee.Reducer.min())
    return image.updateMask(mask.not()).updateMask(mask2a);
    
  };
    
  function applyScaleFactors(image) {
    var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2).multiply(10000);
    return image.addBands(opticalBands, null, true);
  }
  
  // // Aplica um filtro de média morfológica a cada banda de uma imagem usando um kernel personalizado e em seguida é sobreposta pela imagme original
  // function gapfilling(image){
  //   return image.focal_mean(1, 'square', 'pixels', 30).blend(image)
  // }

  var LT05 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
              .filterBounds(roi)
              .filterDate(SD, ED)
              .map(cloudMask_L57)
              .map(applyScaleFactors)
              .select(bands_5,['B1', 'B2', 'B3', 'B4', 'B5', 'B7']);
  
  var LE07 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
              .filterBounds(roi)
              .filterDate(SD, ED)
              .map(cloudMask_L57)
              .map(applyScaleFactors)
              //.map(gapfilling)
              .select(bands_5, ['B1', 'B2', 'B3', 'B4', 'B5', 'B7']);
              
  var LC08 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
              .filterBounds(roi)
              .filterDate(SD, ED)
              .map(cloudMask_L89)
              .map(applyScaleFactors)
              .select(bands_8, ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7']);
  /* 5. Do spectral unmixing and calculate NDFIW */
  // before disturbance
  var unmixNDFIW5 = function(img) {

      var unmixImg = ee.Image(img).unmix([gv_5, shade_5, wat_5, soil_5, cloud_5], true,true) //change the variable #image for the variable you want run (e.g., imageB05, imageA05) 
                    .rename(['GV', 'Shade', 'Water','Soil','Cloud'])
      var imgWithEndmembers = ee.Image(img).addBands(unmixImg)
        
      var ndfiw = imgWithEndmembers.expression(
        '((GV / (1 - SHADE)) - (SOIL + WATER)) / ((GV / (1 - SHADE)) + (SOIL + WATER))', {
            'GV': unmixImg.select('GV'),
          'SHADE': unmixImg.select('Shade'),
            'WATER': unmixImg.select('Water'),
          'SOIL': unmixImg.select('Soil')
        }).rename('NDFIW')
        
        var gvShade = imgWithEndmembers.expression(
        'GV / (1 - SHADE)', {
        'GV': unmixImg.select('GV'),
        'SHADE': unmixImg.select('Shade'),
    }).rename('GVshade')

      return imgWithEndmembers.addBands([ndfiw,gvShade])
  };
  var unmixNDFIW7 = function(img) {

      var unmixImg = ee.Image(img).unmix([gv_7, shade_7, wat_7, soil_7, cloud_7], true,true) //change the variable #image for the variable you want run (e.g., imageB05, imageA05) 
                    .rename(['GV', 'Shade', 'Water','Soil','Cloud'])
      var imgWithEndmembers = ee.Image(img).addBands(unmixImg)
        
      var ndfiw = imgWithEndmembers.expression(
        '((GV / (1 - SHADE)) - (SOIL + WATER)) / ((GV / (1 - SHADE)) + (SOIL + WATER))', {
            'GV': unmixImg.select('GV'),
          'SHADE': unmixImg.select('Shade'),
            'WATER': unmixImg.select('Water'),
          'SOIL': unmixImg.select('Soil')
        }).rename('NDFIW')

      var gvShade = imgWithEndmembers.expression(
        'GV / (1 - SHADE)', {
        'GV': unmixImg.select('GV'),
        'SHADE': unmixImg.select('Shade'),
    }).rename('GVshade')

      return imgWithEndmembers.addBands([ndfiw,gvShade])
  };
  var unmixNDFIW8 = function(img) {

      var unmixImg = ee.Image(img).unmix([gv_8, shade_8, wat_8, soil_8, cloud_8], true,true) //change the variable #image for the variable you want run (e.g., imageB05, imageA05) 
                    .rename(['GV', 'Shade', 'Water','Soil','Cloud'])
      var imgWithEndmembers = ee.Image(img).addBands(unmixImg)
        
      var ndfiw = imgWithEndmembers.expression(
        '((GV / (1 - SHADE)) - (SOIL + WATER)) / ((GV / (1 - SHADE)) + (SOIL + WATER))', {
            'GV': unmixImg.select('GV'),
          'SHADE': unmixImg.select('Shade'),
            'WATER': unmixImg.select('Water'),
          'SOIL': unmixImg.select('Soil')
        }).rename('NDFIW')

      var gvShade = imgWithEndmembers.expression(
        'GV / (1 - SHADE)', {
        'GV': unmixImg.select('GV'),
        'SHADE': unmixImg.select('Shade'),
    }).rename('GVshade')

      return imgWithEndmembers.addBands([ndfiw,gvShade])
  };

  var endmembersCollection05 = LT05.map(unmixNDFIW5)
  var endmembersCollection07 = LE07.map(unmixNDFIW7)
  var endmembersCollection08 = LC08.map(unmixNDFIW8)


  var endmembersCollection_5_7_8 = endmembersCollection05.merge(endmembersCollection07).merge(endmembersCollection08)

  // -------------------------------- Bloco adicionado mask usando cloud fraction -------------------------
  endmembersCollection_5_7_8 = endmembersCollection_5_7_8.map(function(image) {
      var cfThreshold = ee.Image.constant(0.05)
      var mask = image.select('Cloud').lt(cfThreshold)
      
      return ee.Image(image).updateMask(mask).select(['GV','Shade','Water','Soil','NDFIW', 'GVshade'])
      
    })
      
  var trainNDFI = endmembersCollection_5_7_8.filterDate(start, end)

  var medianEM = trainNDFI.select([
     'GVshade',
     'Shade',
     'Water',
     'Soil'
     ]).median()

  var minNDFI = trainNDFI.select('NDFIW').median() // TODO: This could be low

  var median = medianEM.addBands(minNDFI)
  
  var medianMetaData = ee.Image(median).setMulti({
                'month': start.format('M'),
                'year': ee.Number(start.format('Y')),
                'system:time_start': start.millis(),
                'bands': median.bandNames().length()
                })
  
  return medianMetaData
  
} 

exports.getInputsNoRegion = function(start, end, params) {
  /* Get inputs for no specific region */
  
  var SD = start;
  var ED = end;

  // bandas landsat 5
  var bands_5 = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'];
  // bandas landsat 8
  var bands_8 = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']; 
  
  // Endmembers LT05 e LC08

  var cloud_5 = [16022.125,	7781.473195,	7843.746653,	8623.887614,	6047.753435,	3765.897382]
  var gv_5 = [261.5757123,	577.9617759,	366.6097816,	5291.62735,	2121.13877,	770.4923314]
  var soil_5 = [1036.654025,	1466.611126,	1777.302094,	2795.225033,	4631.590805,	3619.778338]
  var wat_5 = [2090.036621,	2955.367136,	1767.187081,	232.7541696,	35.7418541,	21.8578176]
  var shade_5 = [0, 0, 0, 0, 0, 0];
  
  var cloud_7 = [16022.125,	7781.473195,	7843.746653,	8623.887614,	6047.753435,	3765.897382]
  var gv_7 = [261.5757123,	577.9617759,	366.6097816,	5291.62735,	2121.13877,	770.4923314]
  var soil_7 = [1036.654025,	1466.611126,	1777.302094,	2795.225033,	4631.590805,	3619.778338]
  var wat_7 = [2090.036621,	2955.367136,	1767.187081,	232.7541696,	35.7418541,	21.8578176]
  var shade_7 = [0, 0, 0, 0, 0, 0];
  
  var cloud_8 = [ 6375.76141,	6693.028614,	6780.477957,	6968.257616,	7492.911619,	5969.836111,	4662.671087];
  var gv_8 = [281.2302661,	281.1218856,	645.1652039,	269.2813137,	6280.259655,	2007.282559,	735.9335787];
  var soil_8 = [ 1633.97976,	1923.408204,	2963.690719,	3939.457246,	5013.318084,	6035.979192,	4142.394072];
  var wat_8 = [1118.224128,	1611.110856,	2307.231735,	1311.763986,	113.4761094,	181.8146989,	145.5534865];
  var shade_8 = [0, 0, 0, 0, 0, 0, 0];
  
  
  // filtro de nuvem cloudMaskL457
var cloudMask_L57 = function(image) {
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
  var dilated_cloud = (1 << 1);
  var cloud = (1 << 3);
  var cloud_shadow = (1 << 4);
  var snow = (1 << 5);
  var water = (1 << 7);
  var cloud_conf = (1 << 9);
  var cloud_shadow_conf = (1 << 11);
  var snow_ice_conf = (1 << 12);
  var cirrus_conf = (1 << 15);
  
  // Get the pixel QA band.
  var qa = image.select('QA_PIXEL').rename('clouds');
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloud)
  .and(qa.bitwiseAnd(cloud_conf))
  .or((qa.bitwiseAnd(cloud_shadow))
  .and(qa.bitwiseAnd(cloud_shadow_conf)))
  .or((qa.bitwiseAnd(snow))
  .and(qa.bitwiseAnd(snow_ice_conf)))
  .or(qa.bitwiseAnd(cirrus_conf))
  .or(qa.bitwiseAnd(dilated_cloud))
  .or(qa.bitwiseAnd(water))
  var mask2a = image.mask().reduce(ee.Reducer.min())
    return image.updateMask(mask.not()).updateMask(mask2a);
     
};

// filtro de nuvem cloudMaskL457
var cloudMask_L89 = function(image) {
    // Bits 3 and 5 are cloud shadow and cloud, respectively.
    var dilated_cloud = (1 << 1);
    var cirrus = (1 << 2);
    var cloud = (1 << 3);
    var cloud_shadow = (1 << 4);
    var snow = (1 << 5);
    var water = (1 << 7);
    var cloud_conf = (1 << 9);
    var cloud_shadow_conf = (1 << 11);
    var snow_ice_conf  = (1 << 12);
    var cirrus_conf = (1 << 15);
    
    // Get the pixel QA band.
    var qa = image.select('QA_PIXEL');
    // Both flags should be set to zero, indicating clear conditions.
    var mask = qa.bitwiseAnd(cloud)
                  .and(qa.bitwiseAnd(cloud_conf))
                  .or((qa.bitwiseAnd(cloud_shadow))
                  .and(qa.bitwiseAnd(cloud_shadow_conf)))
                  .or((qa.bitwiseAnd(snow))
                  .and(qa.bitwiseAnd(snow_ice_conf)))
                  .or((qa.bitwiseAnd(cirrus))
                  .and(qa.bitwiseAnd(cirrus_conf)))
                  .or(qa.bitwiseAnd(dilated_cloud))
                  .or(qa.bitwiseAnd(water))
    var mask2a = image.mask().reduce(ee.Reducer.min())
    return image.updateMask(mask.not()).updateMask(mask2a);
    
  };
    
  function applyScaleFactors(image) {
    var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2).multiply(10000);
    return image.addBands(opticalBands, null, true);
  }
  // // Aplica um filtro de média morfológica a cada banda de uma imagem usando um kernel personalizado e em seguida é sobreposta pela imagme original
  // function gapfilling(image){
  //   return image.focal_mean(1, 'square', 'pixels', 30).blend(image)
  // }

  var LT05 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
              .filterDate(SD, ED)
              .map(cloudMask_L57)
              .map(applyScaleFactors)
              .select(bands_5,['B1', 'B2', 'B3', 'B4', 'B5', 'B7']);
  var LE07 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
              .filterDate(SD, ED)
              .map(cloudMask_L57)
              .map(applyScaleFactors)
              //.map(gapfilling)
              .select(bands_5,['B1', 'B2', 'B3', 'B4', 'B5', 'B7']);
  var LC08 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
              .filterDate(SD, ED)
              .map(cloudMask_L89)
              .map(applyScaleFactors)
              .select(bands_8, ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7']);
  
  /* 5. Do spectral unmixing and calculate NDFIW */
  // before disturbance
  var unmixNDFIW5 = function(img) {

      var unmixImg = ee.Image(img).unmix([gv_5, shade_5, wat_5, soil_5, cloud_5], true,true) //change the variable #image for the variable you want run (e.g., imageB05, imageA05) 
                    .rename(['GV', 'Shade', 'Water','Soil','Cloud'])
      var imgWithEndmembers = ee.Image(img).addBands(unmixImg)
        
      var ndfiw = imgWithEndmembers.expression(
        '((GV / (1 - SHADE)) - (SOIL + WATER)) / ((GV / (1 - SHADE)) + (SOIL + WATER))', {
            'GV': unmixImg.select('GV'),
          'SHADE': unmixImg.select('Shade'),
            'WATER': unmixImg.select('Water'),
          'SOIL': unmixImg.select('Soil')
        }).rename('NDFIW')
        
        var gvShade = imgWithEndmembers.expression(
        'GV / (1 - SHADE)', {
        'GV': unmixImg.select('GV'),
        'SHADE': unmixImg.select('Shade'),
    }).rename('GVshade')

      return imgWithEndmembers.addBands([ndfiw,gvShade])
  };
  
  var unmixNDFIW7 = function(img) {

      var unmixImg = ee.Image(img).unmix([gv_7, shade_7, wat_7, soil_7, cloud_7], true,true) //change the variable #image for the variable you want run (e.g., imageB05, imageA05) 
                    .rename(['GV', 'Shade', 'Water','Soil','Cloud'])
      var imgWithEndmembers = ee.Image(img).addBands(unmixImg)
        
      var ndfiw = imgWithEndmembers.expression(
        '((GV / (1 - SHADE)) - (SOIL + WATER)) / ((GV / (1 - SHADE)) + (SOIL + WATER))', {
            'GV': unmixImg.select('GV'),
          'SHADE': unmixImg.select('Shade'),
            'WATER': unmixImg.select('Water'),
          'SOIL': unmixImg.select('Soil')
        }).rename('NDFIW')

      var gvShade = imgWithEndmembers.expression(
        'GV / (1 - SHADE)', {
        'GV': unmixImg.select('GV'),
        'SHADE': unmixImg.select('Shade'),
    }).rename('GVshade')

      return imgWithEndmembers.addBands([ndfiw,gvShade])
  };
  
  var unmixNDFIW8 = function(img) {

      var unmixImg = ee.Image(img).unmix([gv_8, shade_8, wat_8, soil_8, cloud_8], true,true) //change the variable #image for the variable you want run (e.g., imageB05, imageA05) 
                    .rename(['GV', 'Shade', 'Water','Soil','Cloud'])
      var imgWithEndmembers = ee.Image(img).addBands(unmixImg)
        
      var ndfiw = imgWithEndmembers.expression(
        '((GV / (1 - SHADE)) - (SOIL + WATER)) / ((GV / (1 - SHADE)) + (SOIL + WATER))', {
            'GV': unmixImg.select('GV'),
          'SHADE': unmixImg.select('Shade'),
            'WATER': unmixImg.select('Water'),
          'SOIL': unmixImg.select('Soil')
        }).rename('NDFIW')

      var gvShade = imgWithEndmembers.expression(
        'GV / (1 - SHADE)', {
        'GV': unmixImg.select('GV'),
        'SHADE': unmixImg.select('Shade'),
    }).rename('GVshade')

      return imgWithEndmembers.addBands([ndfiw,gvShade])
  };

      
  var endmembersCollection05 = LT05.map(unmixNDFIW5)
  var endmembersCollection07 = LE07.map(unmixNDFIW7)
  var endmembersCollection08 = LC08.map(unmixNDFIW8)


  var endmembersCollection_5_7_8 = endmembersCollection05.merge(endmembersCollection07).merge(endmembersCollection08)

  // -------------------------------- Bloco adicionado mask usando cloud fraction -------------------------
  endmembersCollection_5_7_8 = endmembersCollection_5_7_8.map(function(image) {
      var cfThreshold = ee.Image.constant(0.05)
      var mask = image.select('Cloud').lt(cfThreshold)
      
      return ee.Image(image).updateMask(mask).select(['GV','Shade','Water','Soil','NDFIW', 'GVshade'])
      
    })
      
  var trainNDFI = endmembersCollection_5_7_8.filterDate(start, end)

  var medianEM = trainNDFI.select([
     'GVshade',
     'Shade',
     'Water',
     'Soil'
     ]).median()

  var minNDFI = trainNDFI.select('NDFIW').median()

  var median = medianEM.addBands(minNDFI)
  
  var medianMetaData = ee.Image(median).setMulti({
                'month': start.format('M'),
                'year': start.format('Y'),
                'system:time_start': start.millis(),
                'bands': median.bandNames().length()
                })
  return medianMetaData
  
} 

exports.getMonthlyComposites = function(list) {
    /* Get monthly composites */
    var date = ee.Date(ee.List(list).get(0))
    var region = ee.FeatureCollection(ee.List(list).get(1))
    var params = ee.Dictionary(ee.List(list).get(2))
    var collection = ee.ImageCollection([
                      exports.getInputsRegion(date, date.advance(1, 'month'), region, params),
                      exports.getInputsRegion(date.advance(1, 'month'), date.advance(2, 'month'), region, params),
                      exports.getInputsRegion(date.advance(2, 'month'), date.advance(3, 'month'), region, params),
                      exports.getInputsRegion(date.advance(3, 'month'), date.advance(4, 'month'), region, params),
                      exports.getInputsRegion(date.advance(4, 'month'), date.advance(5, 'month'), region, params),
                      exports.getInputsRegion(date.advance(5, 'month'), date.advance(6, 'month'), region, params),
                      exports.getInputsRegion(date.advance(6, 'month'), date.advance(7, 'month'), region, params),
                      exports.getInputsRegion(date.advance(7, 'month'), date.advance(8, 'month'), region, params),
                      exports.getInputsRegion(date.advance(8, 'month'), date.advance(9, 'month'), region, params),
                      exports.getInputsRegion(date.advance(9, 'month'), date.advance(10, 'month'), region, params),
                      exports.getInputsRegion(date.advance(10, 'month'), date.advance(11, 'month'), region, params),
                      exports.getInputsRegion(date.advance(11, 'month'), date.advance(12, 'month'), region, params),
                      ])
    return collection 
    
}

exports.concatenateCollections = function(yearlyCol, fullCol){
  /* Utility function for merging collections */
  return ee.ImageCollection(fullCol).merge(ee.ImageCollection(yearlyCol))
}

exports.getMonthlyCompositesTraining = function(list) {
    // Get monthly composites for the training data locations
    var date = ee.Date(ee.List(list).get(0))
    var params = ee.Dictionary(ee.List(list).get(1))
    
    var collection = ee.ImageCollection([
                      exports.getInputsNoRegion(date, date.advance(1, 'month'), params),
                      exports.getInputsNoRegion(date.advance(1, 'month'), date.advance(2, 'month'), params),
                      exports.getInputsNoRegion(date.advance(2, 'month'), date.advance(3, 'month'), params),
                      exports.getInputsNoRegion(date.advance(3, 'month'), date.advance(4, 'month'), params),
                      exports.getInputsNoRegion(date.advance(4, 'month'), date.advance(5, 'month'), params),
                      exports.getInputsNoRegion(date.advance(5, 'month'), date.advance(6, 'month'), params),
                      exports.getInputsNoRegion(date.advance(6, 'month'), date.advance(7, 'month'), params),
                      exports.getInputsNoRegion(date.advance(7, 'month'), date.advance(8, 'month'), params),
                      exports.getInputsNoRegion(date.advance(8, 'month'), date.advance(9, 'month'), params),
                      exports.getInputsNoRegion(date.advance(9, 'month'), date.advance(10, 'month'), params),
                      exports.getInputsNoRegion(date.advance(10, 'month'), date.advance(11, 'month'), params),
                      exports.getInputsNoRegion(date.advance(11, 'month'), date.advance(12, 'month'), params),
                      ])
    return collection 
    
}

exports.getYears = function(trainBeginning, monitorEnd, region, params){
  /* Utility function for getting a list with study years, region, and parameters */
  var years = ee.List.sequence(trainBeginning, monitorEnd).map(function(int) {
                return ee.Date(ee.String(int).slice(0,4))
  })

  var yearsRegion = years.map(function(year) {
                     return ee.List([year, region, params])
  })
  
  return yearsRegion
}



