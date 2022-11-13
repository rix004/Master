% Get plots
for i = 11:12
    if i == 1
        
        Plots('TPFAConvergenceVoronoi','./Filer/L2errorData')

    elseif i == 2

        Plots('TotalSystemTest','./Filer/L2errorDataSystemTest1');

    elseif i == 3

        Plots('ImpactFieldAndDeterministicTree','TestDETtree')

    elseif i == 4

        Plots('ImpactFieldAndUnstructuredTree','TestDLAtree')

    elseif i == 5

        Plots('ImpactFieldAndUnstructuredTree','TestRRTtree')

    elseif i == 6

        Plots('ErrorPlot_KD_deltam','./Filer/VaryDeltam_and_KD')
        
    elseif i == 7

        Plots('ErrorPlot_alfaR_KD','./Filer/VaryAlphaR_and_KD')

    elseif i == 8

        Plots('ErrorPlot_M_deltaM','./Filer/VaryDeltam_and_M')

    elseif i == 9

        Plots('ErrorPlot_alfaR_M','./Filer/VaryAlphaR_and_M')

    elseif i == 10

        Plots('ErrorPlot_realradii','./Filer/KTvsR_Test_TrueRadii')

    elseif i == 11

        Plots('TrueVScodeKT','./Filer/AnalyticVScompKT1')

    elseif i == 12

        Plots('KTradiusrateFIG','./Filer/KTvsDeltaRData')

    elseif i == 13

        Plots('KTvsR')

    elseif i == 14

        Plots('PressureAndDeviation')

    elseif i == 15

        Plots('KTvsRDLAfig')

    elseif i == 16

        Plots('aVSparticlesRRT')

    elseif i == 17

        Plots('ErrorPlotDifferentValuesDLA','./Filer/DLAerror_KD_and_Nparticles')

    elseif i == 18

        Plots('PressurePlotNetworkDarcy','./Filer/PressurePlotDataDET')

    elseif i == 19

        Plots('PeacemanCorrection',0)

    end
    
end