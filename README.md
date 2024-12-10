# Mangrove Degradation Monitoring

This script detects mangrove degradation using Google Earth Engine (GEE). It analyzes land cover changes from different classes such as mangroves, water, soil, and salt marsh. The purpose is to monitor and detect changes in the mangrove ecosystem over time.

## How to Run

To run this script, follow the steps below:

### Step 1: Set Up Your Google Earth Engine Repository

1. **Create a GEE Repository**:
   - Go to [Google Earth Engine](https://code.earthengine.google.com/).
   - Click on **"New Repository"**.
   - Name your repository (e.g., `Mangrove_Degradation`).

2. **Add Your Scripts**:
   - Upload all the necessary `.js` files (scripts), present in this directory, into your newly created repository.

   - Make sure your `changeDetection.js`, `changeUtils.js`, `classutils.js`, `dataUtils.js`, `miscUtils.js`, `regressionUtils.js`, and `main.js` files are included in the repository.

3. **Import Required Functions**:
   - The script uses external functions from other JavaScript files. You need to adjust the `require` paths for each `.js` file so that they point to the correct location of the required function files.
   - Example for changing `require` paths:

     ```javascript
     var codedUtils = require('users/your_username/your_repository:changeDetection');
     var dataUtils = require('users/your_username/your_repository:dataUtils');
     ```

     Ensure the paths reflect the structure of your GEE repository.

4. **Auto-import from the GEE Repository**:
   - If you want an easy way to load the necessary scripts, you can import the entire repository using the link: 
     [Import Mangrove Degradation Scripts](https://code.earthengine.google.com/?accept_repo=users/igor_cnpy/corescam).

## Step 2: Executing the Monitoring Algorithm `main.js`
### Step 2.1: Defining Parameters

You may adjust the parameters in the script to suit your analysis needs. These parameters control various aspects of the change detection model:

```javascript
var params = ee.Dictionary({
    'cfThreshold': .05,  // Minimum threshold to remove clouds based on cloud fraction
    'consec': 3,  // Consecutive observations beyond change threshold to trigger a change
    'thresh': 3,  // Change threshold defined as an observation's residual normalized by the training models RMSE
    'start': 2000,  // Start year of analysis
    'end': 2022,  // End year of analysis
    'trainDataEnd': 1999,  // End year for the training period
    'trainDataStart': 1997,  // Start year for the training period
    'trainLength': 3,  // Number of years in the training period
    'forestLabel': 1,  // Label for mangrove class
    'window': 2,  // Window size for temporal smoothing
    'minYears': 3,  // Minimum years between disturbances
    'minRMSE': .015,  // Minimum RMSE for generating a valid change score
    'numChanges': 6,  // Maximum number of changes to keep upon export
    'minObs': 5  // Minimum number of observations required for a valid change
});
```

### Step 2.2: Defining the Samples for Training

You have two options for providing the training data:

#### Option 1: Using Points (FeatureCollection)
You can provide the training samples as a `FeatureCollection` in Google Earth Engine. It is essential that each sample is distinct by class and has a unique `id`. The `id` is important because it is used throughout the analysis, especially in the final output map.

Ensure that each class has its own ID to differentiate it properly. For instance, mangroves might have an ID of `1`, water could be `2`, and so on.

#### Option 2: Drawing Samples in Google Earth Engine
If you do not have predefined samples, you can manually draw them in Google Earth Engine using the geometry tools. You can create polygons or points for each land cover class by using the Earth Engine interface.

#### Important Notes:
- **Defining the ID**: It is crucial to define a unique `id` for each class (mangrove, water, soil, salt marsh, etc.). These IDs are used throughout the analysis, including in the final output map.
- **Class IDs in the Script**: In your training data, the IDs for each class should be specified. For example, you may define mangroves with `mangrove: 1`, water with `water: 2`, and so on.

```javascript
// Merge all sample datasets into the training data collection
var trainingData = mangrove.merge(soil).merge(water).merge(saltMarsh);
```

By following these guidelines, you'll ensure that the model can properly distinguish between different land cover classes, and the IDs will be correctly associated throughout the analysis process.

