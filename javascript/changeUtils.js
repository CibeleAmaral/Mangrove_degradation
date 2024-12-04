/* Utilities for main change detection function in CODED */

exports.monitorFunction = function(image, statusList) {
  /* Function to handle iterating over images to calculate consecutive
     observations that look like change */
     
  // Current image in iteration
  var imageUnmasked = ee.Image(image).select('normalizedRes').unmask()

  // Last image in status list
  var previousFull = ee.Image(ee.List(statusList).get(-1))

  // Change threshold on the normalized residuals
  var thresh = previousFull.metadata('thresh') 
  
  // Consecutive observations to flag a change
  var consec = previousFull.metadata('consec')

  // Consecutive observation band from last observation
  var previous = previousFull.select('consecutive'); 
  
  // Change flag from previous iteration
  var changeFlag = previousFull.select('changeFlag');
  
  // Change magnitude from previous iteration
  var previousMagnitude = previousFull.select('magnitude');
  
  // Cloud mask
  var unmaskedObservation = image.select('normalizedRes').mask().eq(1)

  // Current observation beyond threshold and not masked
  var beyondThreshold = imageUnmasked.abs().gte(thresh)
                          .and(unmaskedObservation.eq(1))
                          .and(changeFlag.eq(0))
  
  // Multiplier: 1 for cloud or exceeding threshold, 0 for under threshold
  var multiplierCoef = unmaskedObservation.eq(0).or(beyondThreshold.eq(1))
  
  // Addition factor: 0 for cloud or under threshold, 1 for exceeding
  var addCoef = beyondThreshold.eq(1)
  
  // Add to previous if beyond threshold and not masked
  var consecutive = previous.add(beyondThreshold)
                      .multiply(multiplierCoef)

  // Record magnitude if it is changing 
  var magnitudeNew = imageUnmasked
                      .multiply(beyondThreshold)
                      .rename('magnitude')

  // Flag if consecutive observations exceeds consec parameter
  var changeFlagNew = consecutive.eq(consec).or(changeFlag.eq(1)).rename('changeFlag')
  
  // Get date for every time consecutive = 1. Can filter later. 
  var changeDateNew = consecutive.eq(1)
                        .multiply(image.metadata('system:time_start'))
                        .divide(ee.Image(315576e5)).rename('changeDate') // years since 1970
  
  // Image to add at the end of the status list                                   
  var toAdd = consecutive
              .rename('consecutive')
              .setMulti({'system:time_start': image.get('system:time_start'),
                         'consec': ee.Number(previousFull.get('consec')).short(),
                         'thresh': ee.Number(previousFull.get('thresh')).short()
                        }
                       );
              
  return ee.List(statusList)
           .add((ee.Image(toAdd)
               .addBands([changeFlagNew, magnitudeNew, changeDateNew])
               .cast({'consecutive': 'short',
                      'changeFlag': 'short',
                      'magnitude': 'float',
                      'changeDate': 'float',
                     }
                    )
                )
               )

}

