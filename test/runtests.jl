
using Test

@testset "ch_coordinates" begin

#testing the parameters according to page 19 of 20 of this PDF  https://www.swisstopo.admin.ch/content/swisstopo-internet/de/online/calculation-services/_jcr_content/contentPar/tabs/items/dokumente_und_publik/tabPar/downloadlist/downloadItems/8_1467103085694.download/refsys_d.pdf

#LV03, LN02
have=[602030.680 191775.030 897.915;
617306.300 268507.300 456.064;
776668.105 265372.681 1042.624;
497313.292 145625.438 1207.434;
722758.810 87649.670 1636.600]

#altitudes from step 3
altitudes2=[897.361;
457.138;
1043.616;
1206.367;
1634.472]

have_mod=hcat(have[:,1:2],altitudes2)

want=[    7 27 54.983506 46 52 37.540562 947.149;
7 40 6.983077 47 34 1.385301 504.935;
9 47 3.697723 47 30 55.172797 1089.372;
6 6 7.326361 46 27 14.690021 1258.274;
9 1 16.389053 45 55 45.438020 1685.027]

for i=1:size(want,1)
    thisDatum=have_mod[i,:]
    calculated_result=LV03toETRS89(thisDatum...)

    #check against expected values
    deltah=calculated_result[3]-want[i,7]
    deltacoords=vcat(compound(calculated_result[1]).-want[i,4:6],compound(calculated_result[2]).-want[i,1:3])
    deltas=vcat(deltah,deltacoords)    
    @show deltas
    @test maximum(abs.(deltas))<0.05
    #1 second is roughly 30 meters (for Switzerland)
    #0.05 is roughly 1.5 meters
end

end

