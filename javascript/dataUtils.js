
exports.addDefIter = function(yearImage, statusList) {
  // Iterate through years to:
  // A. find conversion that was not labeled as change but a land cover occured
  // B. Find first year of sustained regrowth.
  
  var lcImage = ee.Image(yearImage)
  var statList = ee.List(statusList)
  var lastImage = ee.Image(statList.get(-1))
  var lastYear = ee.Number(lastImage.get('year'))
  var newYear = ee.Number(lastYear.add(1))
  
  // Make a mask where there was no change recorded in current, previous, or next year
  var changeImage = lastImage.select('dist.*')
  var postImage = lastImage.select('post.*')

  var noChange1 = ee.Image(changeImage.eq(ee.Image(newYear))
                  .or(changeImage.eq(ee.Image(newYear.add(1))))
                  .or(changeImage.eq(ee.Image(newYear.add(2))))
                  .or(changeImage.eq(ee.Image(newYear.subtract(1))))
                  .or(changeImage.eq(ee.Image(newYear.subtract(2))))
                  .reduce(ee.Reducer.sum()))
                  .eq(0)
  
  var noChange = noChange1.eq(1).and(
                  postImage.reduce(ee.Reducer.max()).lt(2)) // No other deforestation
  
  var preChange = lastImage.select('change2.*')
  var preReg = lastImage.select('regrowth2.*')        
  
  // Previously added extra change
  var preAdd = preChange.reduce(ee.Reducer.sum()).gt(0)
  
  // Consecutive water years
  var lastImageWater = lastImage.select('waterConsec')
  var totalWater = lastImage.select('waterTotal')

  // Consecutive forest years
  var lastImageForest = lastImage.select('forestConsec')
  var totalForest = lastImage.select('forestTotal')
  
  // Consecutive nonForest (or water) years
  var lastImageNF = lastImage.select('nfConsec')
  var totalNF = lastImage.select('nfTotal')

  // Flag if no more change is allowed
  // var lastNoChange = lastImage.select('noChange')
  
  var waterValue = ee.Image(2)
  var forestValue = ee.Image(1)
  var zeroImage = ee.Image(0)
  var nnImage = ee.Image(99)
  
  var waterData = lcImage.eq(waterValue)
  var waterConsec = lastImageWater.add(waterData).multiply(waterData).rename('waterConsec')
  var waterTotal = totalWater.add(waterData).rename('waterTotal')
  
  var forestData = lcImage.eq(forestValue)
  var forestConsec = lastImageForest.add(forestData).multiply(forestData).rename('forestConsec')
  var forestTotal = totalForest.add(forestData).rename('forestTotal')

  var nfData = lcImage.neq(forestValue).and(lcImage.neq(zeroImage)).and(lcImage.neq(nnImage)).and(lcImage.neq(waterValue))
  var nfConsec = lastImageNF.add(nfData).multiply(nfData).rename('nfConsec')
  var nfTotal = totalNF.add(nfData).rename('nfTotal')

  var changeBandName = ee.String('change2_').cat(ee.String(newYear))
  var refBandName = ee.String('regrowth2_').cat(ee.String(newYear))

  // Where there was forest, and now classified as non-forest
  var deforestation = forestTotal.gte(ee.Image(4))
                        .and(nfConsec.gte(ee.Image(3)))
                        .multiply(ee.Image(newYear))
                        .multiply(ee.Image(noChange))
                        .multiply(preAdd.eq(0))
                        .rename(changeBandName)

  var reforestation = nfTotal.gte(ee.Image(5))
                        .and(forestTotal.eq(ee.Image(5)))
                        .multiply(ee.Image(newYear))
                        .multiply(ee.Image(noChange))
                        .rename(refBandName)

  var nextImage = changeImage.addBands([nfTotal, nfConsec, forestTotal, forestConsec,
                                       waterTotal, waterConsec, deforestation, reforestation,
                                       preChange, preReg, postImage])
                                       .setMulti({'year': newYear})

  return statList.add(nextImage)
}


                                       


exports.paintFeature = function(feature, width) {
  /* Utility function for converting a feature to an outline */
  var empty = ee.Image().byte();
  var outline = empty.paint({
    featureCollection: ee.FeatureCollection(feature),
    color: 1,
    width: width
  });
  return outline

}

exports.makeBands = function(start, end, prefix) {
  // Make band names for output
  var bandSeq = ee.List.sequence(start, end)

  var bandList = bandSeq.map(function(i) {
    return ee.String(prefix).cat(ee.String(i).slice(0,4))
  })
  
  return bandList
    
}

exports.makeImage = function(arrayImage, column, bandPrefix, start, end) {
  // Turn data array into an image
  var bandList = exports.makeBands(start, end, bandPrefix)
  return arrayImage.arraySlice(1, ee.Number(column), ee.Number(column).add(1)) // (1, 0, 1)
                    .arrayProject([0])
                    .arrayFlatten([bandList])
}

exports.convertCollection = function(iCol) {
  // Convert collection to image with each image in collection becoming a band
  var iterFunc = function(newImage, imageFull) {
    var imageNum = ee.String(ee.Image(imageFull).bandNames().size())
    var bandNames = ee.List(ee.Image(imageFull).bandNames())
                      .add(ee.String('distYear').cat(imageNum))
                      //.add('ah')
    return ee.Image(imageFull)
           .addBands(ee.Image(newImage)).rename(bandNames)

  }
  var startingImage = ee.Image(0)
  var outImage = iCol.iterate(iterFunc, startingImage)
  return ee.Image(outImage)
          .select(
            ee.Image(outImage).bandNames().remove('constant'))
}

exports.convertList = function(iCol) {
  // Convert list to image with each image in collection becoming a band
  var iterFunc = function(newImage, imageFull) {

    return ee.Image(imageFull)
           .addBands(ee.Image(newImage))

  }
  var startingImage = ee.Image(0)
  var outImage = iCol.iterate(iterFunc, startingImage)
  return ee.Image(outImage)
          .select(
            ee.Image(outImage).bandNames().remove('constant'))
}

exports.reduceBands = function(dataImage, params) {
  // Reduce bands to just those containing change. 
  // The 1st argument of the bandNames parameter must be 
  // the prefix that contains the disturbance flag band and
  // the disturbance flag band must be the first set of bands
  // in the image.
  
  var numChanges = ee.Number(params.get('numChanges'))
  var distBands  = exports.makeBands(params.get('start'), params.get('end'), 'dist_')
  var distData = dataImage.select(distBands)
  
  var magBands  = exports.makeBands(params.get('start'), params.get('end'), 'mag_')
  var magData = dataImage.select(magBands)

  var postBands  = exports.makeBands(params.get('start'), params.get('end'), 'post_')
  var postData = dataImage.select(postBands)
  
  var difBands  = exports.makeBands(params.get('start'), params.get('end'), 'dif_')
  var difData = dataImage.select(difBands)

  // Map over years to turn change flag (1) into year of change

  var numYears = ee.Number(params.get('end'))
                .subtract(ee.Number(params.get('start')))
  var numList = ee.List.sequence(0, numYears)
  var yearList = ee.List.sequence(params.get('start'), params.get('end'))

  var changeYears = ee.ImageCollection(numList.map(function(i) {
                      var year = yearList.get(i)
                      var band = distData.select([i])
                      return ee.Image(ee.Number(year))
                              .updateMask(band.eq(1))
                              .rename('year')
                              .toShort()
                    }))
                    
  var changeYearImage = exports.convertCollection(changeYears)
  
  // Need to convert to array to perform accumulate function
  var changeYearArray = distData.unmask().toArray().toArray(1)

  var accum = changeYearArray.arrayAccum(0)
                .arraySlice(1, 0, 1) // (1, 0, 1)
                .arrayProject([0])
                .arrayFlatten([distBands])

  var numList = ee.List.sequence(1, numChanges)

  // Map over every change that is desired and retrieve change information
  var reducedBands = numList.map(function(i) {
                        var distNum = ee.String(ee.Number(i).toInt())
                        var changeYear = ee.Image(changeYearImage
                                            .updateMask(accum.eq(ee.Image(ee.Number(i))))
                                            .reduce(ee.Reducer.min()))
                                            .rename(ee.String('dist_').cat(distNum))
                                            
                        var mag = magData.updateMask(changeYearImage.eq(changeYear))
                                            .reduce(ee.Reducer.min())
                                            .rename(ee.String('mag_').cat(distNum))

                        
                        var post = postData.updateMask(changeYearImage.eq(changeYear))
                                            .reduce(ee.Reducer.min())
                                            .rename(ee.String('post_').cat(distNum))
                                            
                        var dif = difData.updateMask(changeYearImage.eq(changeYear))
                                            .reduce(ee.Reducer.min())
                                            .rename(ee.String('dif_').cat(distNum))

                        return changeYear.addBands([mag, post, dif])
  })

  return exports.convertList(reducedBands)
}

