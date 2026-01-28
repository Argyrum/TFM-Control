load_system("Simulink models\CLRFSPFCRealizedSimscape.slx");
Simulink.exportToVersion("CLRFSPFCRealizedSimscape", "Simulink models\exported\CLRFSPFCRealizedSimscape_2022.slx", "R2022a");
close_system

load_system("Simulink models\OLRippleTest.slx");
Simulink.exportToVersion("OLRippleTest", "Simulink models\exported\OLRippleTest_2022.slx", "R2022a");
close_system

load_system("Simulink models\SSRTest.slx");
Simulink.exportToVersion("SSRTest", "Simulink models\exported\SSRTest_2022.slx", "R2022a");
close_system