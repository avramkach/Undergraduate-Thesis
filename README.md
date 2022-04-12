# Undergraduate-Thesis

Circuit Modeling of Lithium-Ion batteries for a real-time battery management system
• Project goal is to make meaning from online electrochemical impedance measurements done on battery modules and make use of those measurements to make accurate voltage and current prediction of the battery pack
• The work I do to model the batteries directly affects how the battery management system cycles the cells with balancing current limits through the battery modules by using the prediction done by state estimation of SOH, SOC, and Voltages
• Implemented linear Kramer Kronig technique in python to validate EIS measurements
• Designed ECM extraction algorithms that are goings towards a publication
• Made an error model that captures the systems data collection with signal processing
• Using Battery Discharge Dirving profile model is Simulink with Tesla Model S Battery Module EIS and current cycle profile data to analyze voltage predictions from ECM extracted parameters

Folders
• ECM_EXT: Designed Extraction Algorithms scripts
• Error_Model: EIS Measurement Model with function and analysis scripts
• Voltage-Profile-Predictions: Resulting Model S voltage predictions with extracted ECM parameters from EIS data
