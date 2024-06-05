# Data flags and QA names and codes

Below you will find two tables that contain the name and numerical values of flags that can be used to screen measurements for different conditions when using Providentia:

* QA flags relate to GHOST performed quality control checks.
* Data flags relate to standardised flags taken from the data provider.

This information can be used in the configuration files, using either the names of the flags you want to user or their codes.

## QA flags (GHOST checks)

| QA flag name | Code | Default |
| ---      |  ---  | --- |
| Missing Measurement | 0 | ✓ |
| Infinite Value | 1 | ✓ |
| Negative Measurement | 2 | ✓ |
| Zero Measurement | 3 |  |
| Not Maximum Data Quality Level | 4 | |
| Preliminary Data | 5 | |
| Invalid Data Provider Flags - GHOST Decreed | 6 | ✓ |
| Invalid Data Provider Flags - Network Decreed | 7 | |
| No Valid Data to Average | 8 | ✓ |
| Duplicate Station | 9 | ✓ |
| Methodology Not Mapped | 10 | |
| Assumed Primary Sampling | 11 | |
| Assumed Sample Preparation | 12 | |
| Assumed Measurement Methodology | 13 | |
| Unknown Primary Sampling Type | 14 | |
| Unknown Primary Sampling Instrument | 15 | |
| Unknown Sample Preparation Type | 16 | |
| Unknown Sample Preparation Technique | 17 | |
| Unknown Measurement Method | 18 | |
| Unknown Measuring Instrument | 19 | |
| Erroneous Primary Sampling | 20 | ✓ |
| Erroneous Sample Preparation | 21 | ✓ |
| Erroneous Measurement Methodology | 22 | ✓ |
| Invalid QA Measurement Methodology | 23 |
| Sample Gas Volume - Network Standard | 30 |
| Sample Gas Volume - Unknown | 31 |
| Unit Conversion - Network Standard Sample Gas Volume Assumption | 32 | |
| Unit Conversion - Educated Guess Sample Gas Volume Assumption | 33 | |
| Station Position Doubt - DEM Decreed | 40 | |
| Station Position Doubt - Manually Decreed | 41 | |
| Local Precipitation | 50 | |
| Local Extreme Weather | 51 | |
| Local Atmospheric Obscuration | 52 | |
| Local Contamination | 53 | |
| Local Exceptional Event | 54 | |
| Non-Integer Local Timezone (relative to UTC) | 60 | |
| Below Documented Lower Limit of Detection | 70 | |
| Below Reported Lower Limit of Detection | 71 | |
| Below Preferential Lower Limit of Detection | 72 | ✓ |
| Above Documented Upper Limit of Detection | 73 | |
| Above Reported Upper Limit of Detection | 74 | |
| Above Preferential Upper Limit of Detection | 75 | ✓ |
| Insufficient Measurement Resolution - Documented | 80 | |
| Insufficient Measurement Resolution - Reported | 81 | |
| Insufficient Measurement Resolution - Preferential | 82 | ✓ |
| Insufficient Measurement Resolution - Empirical | 83 | ✓ |
| Persistent Recurring Values - 5/6 | 90 | |
| Persistent Recurring Values - 9/12 | 91 | |
| Persistent Recurring Values - 16/24 | 92 | |
| Monthly Fractional Unique Values <= 1% | 100 | |
| Monthly Fractional Unique Values <= 5% | 101 | |
| Monthly Fractional Unique Values <= 10% | 102 | |
| Monthly Fractional Unique Values <= 30% | 103 | |
| Monthly Fractional Unique Values <= 50% | 104 | |
| Monthly Fractional Unique Values <= 70% | 105 | |
| Monthly Fractional Unique Values <= 90% | 106 | |
| Data Outlier - Exceeds Scientifically Decreed Lower/Upper Limit | 110 | ✓ |
| Data Outlier - Monthly Median Exceeds Scientifically Decreed Upper Limit | 111 | ✓ |
| Data Outlier - Network Decreed | 112 | ✓ |
| Data Outlier - Manually Decreed | 113 | ✓ |
| Possible Data Outlier - Monthly Adjusted Boxplot | 114 | |
| Probable Data Outlier - Monthly Adjusted Boxplot | 115 | ✓ |
| Monthly Distribution Consistency - Zone 1 | 120 | |
| Monthly Distribution Consistency - Zone 2 | 121 | |
| Monthly Distribution Consistency - Zone 3 | 122 | |
| Monthly Distribution Consistency - Zone 4 | 123 | |
| Monthly Distribution Consistency - Zone 5 | 124 | |
| Monthly Distribution Consistency - Zone 6 | 125 | |
| Monthly Distribution Consistency - Zone 7 | 126 | |
| Monthly Distribution Consistency - Zone 8 | 127 | |
| Monthly Distribution Consistency - Zone 9 | 128 | |
| Monthly Distribution Consistency - Zone 10 | 129 | | |
| Monthly Distribution Consistency - Unclassified | 130 | |
| Systematic Inconsistent Monthly Distributions - 2/3 Months >= Zone 6 | 131 | |
| Systematic Inconsistent Monthly Distributions - 4/6 Months >= Zone 6 | 132 | ✓ |
| Systematic Inconsistent Monthly Distributions - 8/12 Months >= Zone 6 | 133 | ✓ |

## Data flags (given by data provider)

| Data flag name | Code | Default |
| ---      |  ---  | --- |
|Valid Data | 0 |
|Preliminary Data | 1 |
|Missing Data | 2 |
|Invalid Data - Unspecified | 3 |
|Un-Flagged Data | 4 |
|Estimated Data - Unspecified | 10 |
|Estimated Data - Measured Negative Value | 11 |
|Estimated Data - No Value Detected | 12 |
|Estimated Data - Value Below Detection Limit | 13 |
|Estimated Data - Value Above Detection Limit | 14 |
|Estimated Data - Value Substituted from Secondary Monitor | 15 |
|Estimated Data - Multiple Parameters Aggregated | 16 |
|Extreme/Irregular Data - Unspecified | 20 |
|Data Does Not Meet Internal Network Quality Control Criteria | 21 |
|High Variability of Data | 22 |
|Irregular Data Manually Screened and Accepted | 23 |
|Irregular Data Manually Screened and Rejected | 24 |
|Negative Value | 25 |
|No Value Detected | 26 |
|Reconstructed/Recalculated Data | 27 |
|Value Close to Detection Limit | 28 |
|Value Below Acceptable Range | 29 |
|Value Above Acceptable Range | 30 |
|Value Below Detection Limit | 31 |
|Value Above Detection Limit | 32 |
|Measurement Issue - Unspecified | 40 |
|Chemical Issue | 41 |
|Erroneous Sampling Operation | 42 |
|Extreme Internal Instrument Meteorological Conditions | 43 |
|Extreme Ambient Laboratory Meteorological Conditions | 44 |
|Extreme External Meteorological Conditions | 45 |
|Extreme Sample Transport Conditions | 46 |
|Invalid Flow Rate | 47 |
|Human Error | 48 |
|Matrix Effect | 49 |
|Mechanical Issue/Non-Operational Equipment | 50 |
|No Technician | 51 |
|Operational Maintenance Check Issue | 52 |
|Physical Issue With Filter | 53 |
|Power Failure | 54 |
|Sample Diluted for Analysis | 55 |
|Unmeasured Key Meteorological Parameter | 56 |
|Operational Maintenance - Unspecified | 60 |
|Calibration | 61 |
|Accuracy Check | 62 |
|Blank Check | 63 |
|Detection Limits Check | 64 |
|Precision Check | 65 |
|Retention Time Check | 66 |
|Span Check | 67 |
|Zero Check | 68 |
|Instrumental Inspection | 69 |
|Instrumental Repair | 70 |
|Quality Control Audit | 71 |
|Data Formatting/Processing Issue | 80 |
|Corrected Data Formatting/Processing Issue | 81 |
|Aggregation/Representation Issue - Unspecified | 90 |
|Data Window Completeness < 90% | 91 |
|Data Window Completeness < 75% | 92 |
|Data Window Completeness < 66% | 93 |
|Data Window Completeness < 50% | 94 |
|Data Window Completeness < 25% | 95 |
|>= 75% of Measurements in Window Below Detection Limit | 96 |
|>= 50% of Measurements in Window Below Detection Limit | 97 |
|No Significant Weather | 100 |
|Precipitation - Unspecified Intensity | 101 |
|Precipitation - Light | 102 |
|Precipitation - Moderate | 103 |
|Precipitation - Heavy | 104 |
|Drizzle - Unspecified Intensity | 105 |
|Drizzle - Light | 106 |
|Drizzle - Moderate | 107 |
|Drizzle - Heavy | 108 |
|Freezing Drizzle - Unspecified Intensity | 109 |
|Freezing Drizzle - Light | 110 |
|Freezing Drizzle - Moderate | 111 |
|Freezing Drizzle - Heavy | 112 |
|Rain - Unspecified Intensity | 113 |
|Rain - Light | 114 |
|Rain - Moderate | 115 |
|Rain - Heavy | 116 |
|Rain Shower/s - Unspecified Intensity | 117 |
|Rain Shower/s - Light | 118 |
|Rain Shower/s - Moderate | 119 |
|Rain Shower/s - Heavy | 120 |
|Freezing Rain - Unspecified Intensity | 121 |
|Freezing Rain - Light | 122 |
|Freezing Rain - Moderate | 123 |
|Freezing Rain - Heavy | 124 |
|Freezing Rain Shower/s - Unspecified Intensity | 125 |
|Freezing Rain Shower/s - Light | 126 |
|Freezing Rain Shower/s - Moderate | 127 |
|Freezing Rain Shower/s - Heavy | 128 |
|Snow - Unspecified Intensity | 129 |
|Snow - Light | 130 |
|Snow - Moderate | 131 |
|Snow - Heavy | 132 |
|Snow Shower/s - Unspecified Intensity | 133 |
|Snow Shower/s - Light | 134 |
|Snow Shower/s - Moderate | 135 |
|Snow Shower/s - Heavy | 136 |
|Hail - Unspecified Intensity | 137 |
|Hail - Light | 138 |
|Hail - Moderate | 139 |
|Hail - Heavy | 140 |
|Hail Shower/s - Unspecified Intensity | 141 |
|Hail Shower/s - Light | 142 |
|Hail Shower/s - Moderate | 143 |
|Hail Shower/s - Heavy | 144 |
|Ice Pellets - Unspecified Intensity | 145 |
|Ice Pellets - Light | 146 |
|Ice Pellets - Moderate | 147 |
|Ice Pellets - Heavy | 148 |
|Ice Pellets Shower/s - Unspecified Intensity | 149 |
|Ice Pellets Shower/s - Light | 150 |
|Ice Pellets Shower/s - Moderate | 151 |
|Ice Pellets Shower/s - Heavy | 152 |
|Snow Pellets - Unspecified Intensity | 153 |
|Snow Pellets - Light | 154 |
|Snow Pellets - Moderate | 155 |
|Snow Pellets - Heavy | 156 |
|Snow Pellets Shower/s - Unspecified Intensity | 157 |
|Snow Pellets Shower/s - Light | 158 |
|Snow Pellets Shower/s - Moderate | 159 |
|Snow Pellets Shower/s - Heavy | 160 |
|Snow Grains - Unspecified Intensity | 161 |
|Snow Grains - Light | 162 |
|Snow Grains - Moderate | 163 |
|Snow Grains - Heavy | 164 |
|Diamond Dust - Unspecified Intensity | 165 |
|Diamond Dust - Light | 166 |
|Diamond Dust - Moderate | 167 |
|Diamond Dust - Heavy | 168 |
|Glaze | 169 |
|Rime | 170 |
|Thunderstorm | 171 |
|Funnel Cloud/s | 172 |
|Squalls | 173 |
|Tropical Cyclone (Cyclone/Hurricane/Typhoon) | 174 |
|Duststorm | 175 |
|Sandstorm | 176 |
|Dust/Sand Whirls | 177 |
|High Winds | 178 |
|No Atmospheric Obscuration | 180 |
|Atmospheric Obscuration - Unknown | 181 |
|Dust | 182 |
|Blowing Dust | 183 |
|Drifting Dust | 184 |
|Sand | 185 |
|Blowing Sand | 186 |
|Drifting Sand | 187 |
|Blowing Snow | 188 |
|Drifting Snow | 189 |
|Fog | 190 |
|Freezing Fog | 191 |
|Ground Fog | 192 |
|Ice Fog | 193 |
|Haze | 194 |
|Mist | 195 |
|Sea Spray | 196 |
|Smoke | 197 |
|Volcanic Ash | 198 |
|No Local Contamination | 199 |
|Local Contamination - Unspecified | 200 |
|Agricultural Contamination | 201 |
|Bird-Dropping Contamination | 202 |
|Construction Contamination | 203 |
|Industrial Contamination | 204 |
|Insect Contamination | 205 |
|Internal Laboratory/Instrument Contamination | 206 |
|Pollen/Leaf Contamination | 207 |
|Traffic Contamination | 208 |
|Exceptional Event - Unspecified | 210 |
|Seismic Activity | 211 |
|Stratospheric Ozone Intrusion | 212 |
|Volcanic Eruptions | 213 |
|Wildfire | 214 |
|Chemical Spill/Industrial Accident | 220 |
|Cleanup After a Major Disaster | 221 |
|Demolition | 222 |
|Fireworks | 223 |
|Infrequent Large Gathering | 224 |
|Terrorist Act | 225 |
|Visibility Distance Unlimited | 230 |
|Ceiling Height Unlimited | 231 |
