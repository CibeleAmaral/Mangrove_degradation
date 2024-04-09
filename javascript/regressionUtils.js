// REGRESSION UTILITIES

var miscUtils = require('users/igor_cnpy/corescam:miscUtils');

var constructBandNames = function(base, list) {
  /* Function to get a sequence of band names for harmonic terms. */

  return ee.List(list).map(function(i) {
    return ee.String(base).cat(ee.Number(i).int());
  });
};
  

var regressionTrend = function(iCol, band_name, iColMonitor, numHarmonics, minRMSE) {
  /* Calculate regression coefficients using a trend term */
  
    // The number of cycles per year to model and variable name to use
  var harmonics = numHarmonics
  var dependent = band_name
      
  // Make a list of harmonic frequencies to model.
  // These also serve as band name suffixes.
  var harmonicFrequencies = ee.List.sequence(1, harmonics);
      
  // Construct lists of names for the harmonic terms.
  var cosNames = constructBandNames('cos_', harmonicFrequencies);
  var sinNames = constructBandNames('sin_', harmonicFrequencies);
      
  var addHarmonics = function(freqs) {
   return function(image) {
    // Make an image of frequencies.
    var frequencies = ee.Image.constant(freqs);
    // This band should represent time in radians.
    var time = ee.Image(image).select('t');
    // Get the cosine terms.
    var cosines = time.multiply(frequencies).cos().rename(cosNames);
    // Get the sin terms.
    var sines = time.multiply(frequencies).sin().rename(sinNames);
    return image.addBands(cosines).addBands(sines);
   };
  };    
      
  // Independent variables.
  var independents = ee.List(['constant', 't'])
    .cat(cosNames).cat(sinNames);
      
  // Function to add a time band.
  var addDependents = function(image) {
  // Compute time in fractional years since the epoch.
    var years = image.date().difference('1970-01-01', 'year');
    var timeRadians = ee.Image(years.multiply(2 * Math.PI)).rename('t');
    var constant = ee.Image(1);
    return image.addBands(constant).addBands(timeRadians.float());
  };
      
  // Add variables.
  var harmonicLandsat = ee.ImageCollection(iCol)
    .map(addDependents)
    .map(addHarmonics(harmonicFrequencies));

      
  // The output of the regression reduction is a 4x1 array image.
  var harmonicTrend = harmonicLandsat
    .select(independents.add(dependent))
    .reduce(ee.Reducer.linearRegression(independents.length(), 1));
      
  // Turn the array image into a multi-band image of coefficients.
  var harmonicTrendCoefficients = harmonicTrend.select('coefficients')
    .arrayProject([0])
    .arrayFlatten([independents]);
    
  var fittedHarmonicTraining = harmonicLandsat.map(function(image) {
    return image.addBands(
      image.select(independents)
        .multiply(harmonicTrendCoefficients)
        .reduce('sum')
        .rename('fitted'));
      });  
    
  // Do it again for monitoring period
  var harmonicLandsatMonitor = ee.ImageCollection(iColMonitor)
    .map(addDependents)
    .map(addHarmonics(harmonicFrequencies));
    
  var harmonicTrendMonitoring = harmonicLandsatMonitor
    .select(independents.add(dependent))
    .reduce(ee.Reducer.linearRegression(independents.length(), 1));
    
  var harmonicTrendCoefficientsMon = harmonicTrendMonitoring.select('coefficients')
    .arrayProject([0])
    .arrayFlatten([independents]);
  
  // Using coefficients from monitoring period  
  var fittedHarmonicTrainingMon = harmonicLandsatMonitor.map(function(image) {
    return image.addBands(
      image.select(independents)
        .multiply(harmonicTrendCoefficientsMon)
        .reduce('sum')
        .rename('fitted'));
      });    
    
   // Using coefficients from training period   
  var fittedHarmonicMonitor = harmonicLandsatMonitor.map(function(image) {
    return image.addBands(
      image.select(independents)
        .multiply(harmonicTrendCoefficients)
        .reduce('sum')
        .rename('fitted'));
      });
      
  // Compute residuals manually
  var calc_residuals = function(img){
    var resid = (img.select(dependent)).subtract(img.select('fitted')).rename('residuals')
    return img.addBands(resid)
    }
  
  // Training residuals
  var residualsTraining = fittedHarmonicTraining.map(calc_residuals).select(['residuals'])

  // Monitoring residuals with coefficients from training
  var residualsMonitor = fittedHarmonicMonitor.map(calc_residuals).select(['residuals'])

  // Monitoring residuals with coefficients from monitoring for CUSUM test
  var residualsMonitorCusum = fittedHarmonicTrainingMon.map(calc_residuals).select(['residuals'])

  
  var squaredResiduals = residualsMonitor.map(function(image) {
    var squaredRes = image.select('residuals').pow(2)
    return image.addBands(squaredRes).rename(['residuals','squaredRes'])
  })
  
  var squaredResidualsTraining = residualsTraining.map(function(image) {
    var squaredRes = image.select('residuals').pow(2)
    return image.addBands(squaredRes).rename(['residuals','squaredRes'])
  })  
  
  var rmseNDFIW = ee.ImageCollection(squaredResidualsTraining).select('squaredRes').mean().sqrt().rename('rmse')
  
  var rmseAboveMin = rmseNDFIW.gte(minRMSE).multiply(rmseNDFIW)
  var rmseBelowMin = rmseNDFIW.lt(minRMSE).multiply(ee.Image(minRMSE))
  
  var rmseNDFIWNormalze = rmseAboveMin.add(rmseBelowMin).rename('rmse')
  
  var residualsRMSE = squaredResiduals.map(function(image) {
    var normalizedRes = image.select('residuals').divide(rmseNDFIWNormalze.select('rmse')).rename('normalizedRes')
    return image.addBands(normalizedRes)
  })
  //return ee.Image(1)
      
  return [harmonicTrendCoefficients.rename(['coef_constant','coef_trend','coef_sin','coef_cos']),
          fittedHarmonicMonitor,
          residualsRMSE,
          rmseNDFIW,
          residualsMonitorCusum]
  
}

var regressionNoTrend = function(iCol, band_name, iColMonitor, numHarmonics, minRMSE) {
  /* Calculate regression coefficients without a trend term */
  
    // The number of cycles per year to model and variable name to use
  var harmonics = numHarmonics
  var dependent = band_name
      
  // Make a list of harmonic frequencies to model.
  // These also serve as band name suffixes.
  var harmonicFrequencies = ee.List.sequence(1, harmonics);
      
  // Construct lists of names for the harmonic terms.
  var cosNames = constructBandNames('cos_', harmonicFrequencies);
  var sinNames = constructBandNames('sin_', harmonicFrequencies);
      
  var addHarmonics = function(freqs) {
   return function(image) {
    // Make an image of frequencies.
    var frequencies = ee.Image.constant(freqs);
    // This band should represent time in radians.
    var time = ee.Image(image).select('t');
    // Get the cosine terms.
    var cosines = time.multiply(frequencies).cos().rename(cosNames);
    // Get the sin terms.
    var sines = time.multiply(frequencies).sin().rename(sinNames);
    return image.addBands(cosines).addBands(sines);
   };
  };    
      
  // Independent variables.
  var independents = ee.List(['constant'])
    .cat(cosNames).cat(sinNames);
      
  // Function to add a time band.
  var addDependents = function(image) {
  // Compute time in fractional years since the epoch.
    var years = image.date().difference('1970-01-01', 'year');
    var timeRadians = ee.Image(years.multiply(2 * Math.PI)).rename('t');
    var constant = ee.Image(1);
    return image.addBands(constant).addBands(timeRadians.float());
  };
      
  // Add variables.
  var harmonicLandsat = ee.ImageCollection(iCol)
    .map(addDependents)
    .map(addHarmonics(harmonicFrequencies))
    
  var harmonicLandsatMonitor = ee.ImageCollection(iColMonitor)
    .map(addDependents)
    .map(addHarmonics(harmonicFrequencies))
      
  // The output of the regression reduction is a 4x1 array image.
  var harmonicTrend = harmonicLandsat
    .select(independents.add(dependent))
    .reduce(ee.Reducer.linearRegression(independents.length(), 1));
      
  // Turn the array image into a multi-band image of coefficients.
  var harmonicTrendCoefficients = harmonicTrend.select('coefficients')
    .arrayProject([0])
    .arrayFlatten([independents]);
    
/*    
  var harmonicTrendMon = harmonicLandsatMonitor
    .select(independents.add(dependent))
    .reduce(ee.Reducer.linearRegression(independents.length(), 1));
      
  // Turn the array image into a multi-band image of coefficients.
  var harmonicTrendCoefficientsMon = harmonicTrendMon.select('coefficients')
    .arrayProject([0])
    .arrayFlatten([independents]);
*/

  var fittedHarmonicTraining = harmonicLandsat.map(function(image) {
    return image.addBands(
      image.select(independents)
        .multiply(harmonicTrendCoefficients)
        .reduce('sum')
        .rename('fitted'));
      });  
    
  var fittedHarmonicMonitor = harmonicLandsatMonitor.map(function(image) {
    return image.addBands(
      image.select(independents)
        .multiply(harmonicTrendCoefficients)
        .reduce('sum')
        .rename('fitted'));
      });
      
  // Compute residuals manually
  var calc_residuals = function(img){
    var resid = (img.select(dependent)).subtract(img.select('fitted')).rename('residuals')
    return img.addBands(resid)
    }
  
  var residualsTraining = fittedHarmonicTraining.map(calc_residuals).select(['residuals'])
  var residualsMonitor = fittedHarmonicMonitor.map(calc_residuals).select(['residuals'])
  
  var squaredResiduals = residualsMonitor.map(function(image) {
    var squaredRes = image.select('residuals').pow(2)
    return image.addBands(squaredRes).rename(['residuals','squaredRes'])
  })
  
  var squaredResidualsTraining = residualsTraining.map(function(image) {
    var squaredRes = image.select('residuals').pow(2)
    return image.addBands(squaredRes).rename(['residuals','squaredRes'])
  })  
  
  var rmseNDFIW = ee.ImageCollection(squaredResidualsTraining).select('squaredRes').mean().sqrt().rename('rmse')

  var rmseAboveMin = rmseNDFIW.gte(ee.Image(minRMSE)).multiply(rmseNDFIW)
  var rmseBelowMin = rmseNDFIW.lt(ee.Image(minRMSE)).multiply(ee.Image(minRMSE))
  
  var rmseNDFIWNormalze = rmseAboveMin.add(rmseBelowMin).rename('rmse')
  
  var residualsRMSE = squaredResiduals.map(function(image) {
    var normalizedRes = image.select('residuals').divide(rmseNDFIWNormalze.select('rmse')).rename('normalizedRes')
    return image.addBands(normalizedRes)
  })
  //return ee.Image(1)
      
  return [harmonicTrendCoefficients.rename(['coef_constant_' + band_name,'coef_sin_' + band_name,'coef_cos_'+ band_name]),
          fittedHarmonicMonitor,
          residualsRMSE,
          rmseNDFIW.rename('rmse_' + band_name),
          residualsMonitor]
  
}


exports.doRegression = function(iCol, band_name, iColMonitor, numHarmonics, trend, minRMSE) {
  /* Decide which regression function to use and do it */
  
  var regressionResults = ee.Algorithms.If(ee.Number(trend).eq(1),
                          regressionTrend(iCol, band_name, iColMonitor, numHarmonics, minRMSE),
                          regressionNoTrend(iCol, band_name, iColMonitor, numHarmonics, minRMSE))

  return regressionResults

}

exports.getSquaredResiduals = function(image) {
  /* Calculate squared residuals */
  var squaredRes = image.select('residuals').pow(2)
  return squaredRes
}


var regressionNoTrendTraining = function(iCol, band_name, numHarmonics, minRMSE) {
  
    // The number of cycles per year to model and variable name to use
  var harmonics = numHarmonics
  var dependent = band_name
      
  // Make a list of harmonic frequencies to model.
  // These also serve as band name suffixes.
  var harmonicFrequencies = ee.List.sequence(1, harmonics);
      
  // Construct lists of names for the harmonic terms.
  var cosNames = constructBandNames('cos_', harmonicFrequencies);
  var sinNames = constructBandNames('sin_', harmonicFrequencies);
      
  var addHarmonics = function(freqs) {
   return function(image) {
    // Make an image of frequencies.
    var frequencies = ee.Image.constant(freqs);
    // This band should represent time in radians.
    var time = ee.Image(image).select('t');
    // Get the cosine terms.
    var cosines = time.multiply(frequencies).cos().rename(cosNames);
    // Get the sin terms.
    var sines = time.multiply(frequencies).sin().rename(sinNames);
    return image.addBands(cosines).addBands(sines);
   };
  };    
      
  // Independent variables.
  var independents = ee.List(['constant'])
    .cat(cosNames).cat(sinNames);
      
  // Function to add a time band.
  var addDependents = function(image) {
  // Compute time in fractional years since the epoch.
    var years = image.date().difference('1970-01-01', 'year');
    var timeRadians = ee.Image(years.multiply(2 * Math.PI)).rename('t');
    var constant = ee.Image(1);
    return image.addBands(constant).addBands(timeRadians.float());
  };
      
  // Add variables.
  var harmonicLandsat = ee.ImageCollection(iCol)
    .map(addDependents)
    .map(addHarmonics(harmonicFrequencies))
      
  // The output of the regression reduction is a 4x1 array image.
  var harmonicTrend = harmonicLandsat
    .select(independents.add(dependent))
    .reduce(ee.Reducer.linearRegression(independents.length(), 1));
      
  // Turn the array image into a multi-band image of coefficients.
  var harmonicTrendCoefficients = harmonicTrend.select('coefficients')
    .arrayProject([0])
    .arrayFlatten([independents]);
    
  var fittedHarmonicTraining = harmonicLandsat.map(function(image) {
    return image.addBands(
      image.select(independents)
        .multiply(harmonicTrendCoefficients)
        .reduce('sum')
        .rename('fitted'));
      });  
      
  // Compute residuals manually
  var calc_residuals = function(img){
    var resid = (img.select(dependent)).subtract(img.select('fitted')).rename('residuals')
    return img.addBands(resid)
    }
  
  var residualsTraining = fittedHarmonicTraining.map(calc_residuals).select(['residuals'])
  
  var squaredResidualsTraining = residualsTraining.map(function(image) {
    var squaredRes = image.select('residuals').pow(2)
    return image.addBands(squaredRes).rename(['residuals','squaredRes'])
  })  
  
  var rmseNDFIW = ee.ImageCollection(squaredResidualsTraining).select('squaredRes').mean().sqrt().rename('rmse')

  var rmseAboveMin = rmseNDFIW.gte(ee.Image(minRMSE)).multiply(rmseNDFIW)
  var rmseBelowMin = rmseNDFIW.lt(ee.Image(minRMSE)).multiply(ee.Image(minRMSE))
  
  var rmseNDFIWNormalze = rmseAboveMin.add(rmseBelowMin).rename('rmse')
  
  var residualsRMSE = squaredResidualsTraining.map(function(image) {
    var normalizedRes = image.select('residuals').divide(rmseNDFIWNormalze.select('rmse')).rename('normalizedRes')
    return image.addBands(normalizedRes)
  })
      
  return ee.Image(harmonicTrendCoefficients.addBands(rmseNDFIW).rename(
            ['coef_constant_' + band_name,
            'coef_sin_' + band_name,
            'coef_cos_' + band_name,
            'rmse_' + band_name]))
  
}



exports.doRegressionTraining = function(iCol, band_name, numHarmonics, trend, minRMSE) {
  
  var regressionResults = regressionNoTrendTraining(iCol, band_name, numHarmonics, minRMSE)
  return regressionResults
}


exports.allRegression = function(inputsTraining, inputsMonitoring, numHarmonics, trend, minimumRMSE) {
  
  // Get regression coefficients for all endmember fractions
  
  // 1.1.1 NDFIW
  var minRMSE = ee.Number(minimumRMSE)

  var regressionNDFIW = exports.doRegression(inputsTraining, 'NDFIW', inputsMonitoring, 1, 0, minRMSE)

  var coefficientsImageNDFIW = ee.Image(ee.List(regressionNDFIW).get(0))
  var fittedNDFIW = ee.List(regressionNDFIW).get(1)
  var residualsNDFIW = ee.ImageCollection(ee.List(regressionNDFIW).get(2))
  var rmseNDFIW = ee.Image(ee.List(regressionNDFIW).get(3))
  
  // 1.1.2 GVshade
  
  var regressionGV = exports.doRegression(inputsTraining, 'GVshade', inputsMonitoring, 1, 0, minRMSE)
  var coefficientsImageGV = ee.Image(ee.List(regressionGV).get(0))
  var rmseGV = ee.Image(ee.List(regressionGV).get(3))

  // 1.1.3 Water
  
  var regressionWater = exports.doRegression(inputsTraining, 'Water', inputsMonitoring, 1, 0, minRMSE)
  var coefficientsImageWater = ee.Image(ee.List(regressionWater).get(0))
  var rmseWater = ee.Image(ee.List(regressionWater).get(3))

  // 1.1.4 Soil
  
  var regressionSoil = exports.doRegression(inputsTraining, 'Soil', inputsMonitoring, 1, 0, minRMSE)
  var coefficientsImageSoil = ee.Image(ee.List(regressionSoil).get(0))
  var rmseSoil = ee.Image(ee.List(regressionSoil).get(3))

  // 1.1.5 Shade
  
  var regressionShade = exports.doRegression(inputsTraining, 'Shade', inputsMonitoring, 1, 0, minRMSE)
  var coefficientsImageShade = ee.Image(ee.List(regressionShade).get(0))
  var rmseShade = ee.Image(ee.List(regressionShade).get(3))

  var rmses = rmseNDFIW.addBands([rmseGV, rmseWater, rmseSoil, rmseShade])
  var coefs = coefficientsImageNDFIW.addBands([coefficientsImageGV,coefficientsImageWater,coefficientsImageSoil,coefficientsImageShade])
  var allCoefs = coefs.addBands(rmses)
  
  // 2.1 NDFIW for CUSUM
    
  var regressionNDFIWCu = exports.doRegression(inputsMonitoring, 'NDFIW', inputsMonitoring, 1, 1, minRMSE)
  var residualsNDFIWCu = ee.ImageCollection(ee.List(regressionNDFIWCu).get(4))

  return [residualsNDFIW, allCoefs, rmses, residualsNDFIWCu]

}

// Calculate scaled cumulative sum of residuals. Takes residuals from 
// OLS regression and returns the scaled version, along with a breakpoint
// band. Modified from code from Paulo Arevalo.

exports.squared = function(img){
  /* Get image squared */
  // Modified from Paulo Arevalo's code.
  var imgsquared = img.multiply(img)
                    .setMulti({'system:time_start': img.get('system:time_start'),
                                'year': ee.Number(img.get('year'))})
  return imgsquared
}


exports.calc_scaled_cusum = function(resids, cusum, binaryClassifier, iCol, params, nobs) {
  /* Calculate CUSUM change probability based on OLS Residuals */
  // Scaled cumulative residuals, based on statsmodels (Ploberger, Werner, and Walter Kramer 1992)
    
  var nobssigma2 = resids.map(exports.squared).sum()
  var ddof = ee.Number(4) //DANGER, hardcoded bc code was refactored into functions
  nobssigma2 = nobssigma2.divide(nobs.subtract(ddof)).multiply(nobs)
  
  // Actual calculation of scaled CUSUM. Append breakpoint band
  var scaled_resid = cusum.map(function(img){
    var scaledres = img.divide(nobssigma2.sqrt()).set('system:time_start', img.get('system:time_start')).set('year', ee.Number(img.get('year')))
    var breakpoint = img.metadata('system:time_start').divide(ee.Image(315576e5)).rename('dateBand')
    return scaledres.abs().addBands(breakpoint) 
  })
  
  var maxCusum = scaled_resid.qualityMosaic('cusum').rename(['cusum','cusumDate'])
  
  // Mask before the change and classify coefficients to get binary probability
  
  var afterCusum = ee.List(misc_utils.GetInputsRetrain(iCol, maxCusum.select('cusumDate'), 1, params))
  var afterCusumICol = ee.ImageCollection(afterCusum.get(0))
  var afterCusumReg =  ee.List(exports.allRegression(afterCusumICol, afterCusumICol, 1, 0, params.get('minRMSE')))

  var afterCusumCoefs = ee.Image(afterCusumReg.get(1))
  var afterCusumForestProb = afterCusumCoefs.classify(binaryClassifier)
  
  return maxCusum.addBands(afterCusumForestProb).rename(['cusum','cusumDate','postCuFProb'])
  
}
