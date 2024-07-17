set terminal png small
set output "/Users/test/Desktop/UB/Advanced_Bioinformatics/exercise_02/blast/chrI-x-SRR6130428/chrI-x-SRR6130428_megablast.dnadiff.png"
set ytics ( \
 "7316" 1, \
 "6320" 1173, \
 "*6124" 1469, \
 "*4218" 1721, \
 "5286" 1844, \
 "5552" 2008, \
 "*4874" 2197, \
 "*4490" 2332, \
 "*4038" 2458, \
 "*6670" 2575, \
 "*8176" 2977, \
 "8348" 7189, \
 "5158" 12286, \
 "5760" 12440, \
 "*3956" 12648, \
 "*4860" 12763, \
 "*6436" 12896, \
 "*6964" 13217, \
 "5280" 13826, \
 "*5698" 13989, \
 "7516" 14190, \
 "5878" 15923, \
 "6878" 16144, \
 "*7506" 16656, \
 "6608" 18363, \
 "6850" 18735, \
 "*5234" 19236, \
 "6694" 19396, \
 "*6938" 19812, \
 "5656" 20381, \
 "6798" 20578, \
 "*5768" 21039, \
 "*4940" 21248, \
 "7108" 21386, \
 "*6216" 22198, \
 "*7032" 22472, \
 "6954" 23173, \
 "5960" 23767, \
 "*3857" 23997, \
 "*3920" 24109, \
 "*5650" 24223, \
 "6460" 24419, \
 "4914" 24747, \
 "*5890" 24883, \
 "6494" 25106, \
 "*5446" 25446, \
 "6012" 25626, \
 "5598" 25861, \
 "5292" 26053, \
 "4970" 26218, \
 "*6836" 26358, \
 "*6260" 26838, \
 "8718" 27119, \
 "8646" 34698, \
 "7096" 41687, \
 "*7846" 42484, \
 "*9456" 45311, \
 "*7786" 65334, \
 "*7510" 67884, \
 "*6204" 69614, \
 "7118" 69885, \
 "6990" 70709, \
 "*7996" 71364, \
 "9378" 74856, \
 "*8822" 92461, \
 "8698" 101135, \
 "*9394" 108602, \
 "6014" 126574, \
 "*7334" 126809, \
 "9042" 128036, \
 "*8370" 139157, \
 "*6868" 144390, \
 "*6612" 144897, \
 "5232" 145273, \
 "*7208" 145432, \
 "7332" 146435, \
 "*7364" 147654, \
 "7350" 148952, \
 "*8684" 150224, \
 "*7818" 157555, \
 "*8408" 160281, \
 "4672" 165730, \
 "*3899" 165856, \
 "*4574" 165969, \
 "4730" 166095, \
 "*3859" 166221, \
 "3985" 166333, \
 "4600" 166449, \
 "*4017" 166575, \
 "*4566" 166692, \
 "*4514" 166818, \
 "*3934" 166944, \
 "*4488" 167058, \
 "*3891" 167184, \
 "4025" 167297, \
 "*4166" 167414, \
 "*4686" 167535, \
 "*4710" 167661, \
 "5516" 167787, \
 "*4916" 167973, \
 "4284" 168110, \
 "3991" 168235, \
 "*5410" 168351, \
 "5018" 168526, \
 "5110" 168669, \
 "*4780" 168818, \
 "5368" 168946, \
 "*4550" 169116, \
 "5422" 169242, \
 "*5166" 169418, \
 "4628" 169573, \
 "*5434" 169699, \
 "*3979" 169877, \
 "*3817" 169993, \
 "*4350" 170104, \
 "*4696" 170230, \
 "5064" 170356, \
 "*4736" 170502, \
 "*5826" 170628, \
 "4272" 170843, \
 "*4242" 170967, \
 "4450" 171091, \
 "6158" 171217, \
 "*5144" 171479, \
 "*5408" 171631, \
 "*4552" 171806, \
 "*6660" 171932, \
 "*4398" 172331, \
 "4380" 172457, \
 "5190" 172583, \
 "*5732" 172740, \
 "*5336" 172946, \
 "*4034" 173114, \
 "*3987" 173231, \
 "*4362" 173347, \
 "5920" 173473, \
 "*8456" 173700, \
 "8280" 179470, \
 "*4770" 184261, \
 "*6214" 184388, \
 "8182" 184662, \
 "*7394" 188920, \
 "4958" 190294, \
 "5078" 190433, \
 "*4454" 190580, \
 "6352" 190706, \
 "6358" 191009, \
 "*4982" 191313, \
 "5746" 191453, \
 "5302" 191660, \
 "8082" 191826, \
 "*6398" 195608, \
 "6880" 195920, \
 "5588" 196433, \
 "*6100" 196624, \
 "*6552" 196872, \
 "6638" 197228, \
 "6376" 197613, \
 "*5652" 197920, \
 "*5404" 198117, \
 "7876" 198291, \
 "8688" 201241, \
 "*8030" 208605, \
 "6550" 212225, \
 "6830" 212581, \
 "4114" 213058, \
 "4570" 213178, \
 "*5180" 213304, \
 "6128" 213460, \
 "*5668" 213712, \
 "*5000" 213910, \
 "6024" 214051, \
 "6662" 214288, \
 "6328" 214687, \
 "*5558" 214984, \
 "*6416" 215173, \
 "4978" 215488, \
 "5164" 215628, \
 "*6496" 215783, \
 "*5382" 216123, \
 "6212" 216295, \
 "6450" 216569, \
 "5924" 216895, \
 "*5900" 217122, \
 "6914" 217345, \
 "*4212" 217893, \
 "*5874" 218016, \
 "4108" 218237, \
 "7098" 218356, \
 "5994" 219153, \
 "*5430" 219385, \
 "5770" 219562, \
 "*6074" 219771, \
 "*6802" 220015, \
 "*5742" 220478, \
 "6936" 220684, \
 "6702" 221252, \
 "5236" 221670, \
 "6840" 221830, \
 "*6600" 222318, \
 "*7590" 222687, \
 "6842" 224602, \
 "*5876" 225094, \
 "7522" 225315, \
 "*5270" 227055, \
 "5996" 227217, \
 "*7436" 227450, \
 "*5836" 228937, \
 "*7808" 229153, \
 "*5582" 231839, \
 "*6946" 232030, \
 "7264" 232600, \
 "5394" 233670, \
 "*6160" 233843, \
 "7504" 234106, \
 "4102" 235800, \
 "7366" 235919, \
 "6400" 237226, \
 "*6432" 237538, \
 "*6698" 237858, \
 "*4902" 238275, \
 "*7336" 238411, \
 "*7088" 239642, \
 "5388" 240427, \
 "6692" 240600, \
 "6896" 241014, \
 "6566" 241547, \
 "8858" 241908, \
 "6784" 250964, \
 "7544" 251419, \
 "5392" 253188, \
 "6030" 253361, \
 "4452" 253598, \
 "8830" 253724, \
 "7112" 262467, \
 "8910" 263282, \
 "9334" 272884, \
 "4596" 289031, \
 "5560" 289157, \
 "7228" 289346, \
 "4586" 290360, \
 "7014" 290486, \
 "4003" 291171, \
 "5730" 291287, \
 "4674" 291492, \
 "4568" 291618, \
 "5222" 291744, \
 "5188" 291903, \
 "5084" 292060, \
 "9290" 292207, \
 "5118" 307406, \
 "5182" 307556, \
 "7478" 307713, \
 "7398" 309320, \
 "6354" 310709, \
 "6834" 311012, \
 "6322" 311492, \
 "4296" 311788, \
 "6582" 311913, \
 "6598" 312277, \
 "4390" 312644, \
 "5622" 312770, \
 "7440" 312963, \
 "6700" 314455, \
 "4274" 314873, \
 "5654" 314998, \
 "6736" 315195, \
 "5510" 315631, \
 "4688" 315817, \
 "9644" 315943, \
 "9522" 346481, \
 "4344" 368763, \
 "3770" 368889, \
 "6046" 368998, \
 "6190" 369239, \
 "4122" 369508, \
 "5374" 369628, \
 "6772" 369799, \
 "7368" 370251, \
 "6696" 371564, \
 "6470" 371980, \
 "5120" 372310, \
 "6888" 372460, \
 "6932" 372983, \
 "4328" 373547, \
 "7196" 373673, \
 "6950" 374657, \
 "5142" 375244, \
 "9672" 375396, \
 "5330" 409129, \
 "4646" 409297, \
 "6570" 409423, \
 "9102" 409784, \
 "5148" 421595, \
 "5636" 421748, \
 "5864" 421942, \
 "5600" 422162, \
 "5992" 422354, \
 "6300" 422586, \
 "9550" 422879, \
 "4138" 446685, \
 "4870" 446805, \
 "4062" 446939, \
 "5500" 447057, \
 "7864" 447241, \
 "4256" 450130, \
 "5614" 450254, \
 "4980" 450446, \
 "4608" 450586, \
 "8866" 450712, \
 "4021" 459868, \
 "5348" 459985, \
 "9238" 460154, \
 "7850" 474571, \
 "7142" 477417, \
 "4140" 478309, \
 "7236" 478429, \
 "6486" 479449, \
 "7578" 479786, \
 "4606" 481657, \
 "5102" 481783, \
 "5962" 481932, \
 "8582" 482162, \
 "4808" 488695, \
 "4468" 488825, \
 "6150" 488951, \
 "4252" 489210, \
 "5934" 489334, \
 "6062" 489562, \
 "5272" 489805, \
 "7302" 489968, \
 "4658" 491106, \
 "7550" 491232, \
 "5238" 493011, \
 "6086" 493171, \
 "5454" 493416, \
 "4384" 493597, \
 "6980" 493723, \
 "7156" 494355, \
 "6172" 495271, \
 "4524" 495536, \
 "4248" 495662, \
 "4530" 495786, \
 "4100" 495912, \
 "5010" 496031, \
 "4486" 496173, \
 "5664" 496299, \
 "8006" 496497, \
 "4005" 500035, \
 "6656" 500151, \
 "5762" 500546, \
 "5572" 500754, \
 "5606" 500944, \
 "7248" 501136, \
 "4346" 502180, \
 "8064" 502306, \
 "5184" 506026, \
 "4894" 506183, \
 "5092" 506319, \
 "4538" 506467, \
 "8134" 506593, \
 "4876" 510625, \
 "5300" 510760, \
 "4011" 510925, \
 "7524" 511042, \
 "6838" 512784, \
 "4156" 513270, \
 "5944" 513391, \
 "5230" 513620, \
 "9482" 513779, \
 "6988" 534883, \
 "6614" 535525, \
 "4540" 535901, \
 "6792" 536027, \
 "5846" 536484, \
 "4428" 536702, \
 "5648" 536828, \
 "3781" 537024, \
 "6490" 537133, \
 "4232" 537471, \
 "4578" 537594, \
 "8078" 537720, \
 "6022" 541493, \
 "4664" 541729, \
 "5416" 541855, \
 "4572" 542031, \
 "5670" 542157, \
 "6276" 542355, \
 "8724" 542640, \
 "6796" 550271, \
 "5646" 550731, \
 "5372" 550927, \
 "6968" 551098, \
 "6102" 551711, \
 "4832" 551959, \
 "4324" 552090, \
 "8706" 552216, \
 "5068" 559736, \
 "4998" 559882, \
 "3897" 560023, \
 "5340" 560136, \
 "6640" 560304, \
 "3989" 560689, \
 "5780" 560805, \
 "4838" 561015, \
 "4766" 561147, \
 "3999" 561274, \
 "4358" 561390, \
 "5870" 561516, \
 "6308" 561737, \
 "4342" 562031, \
 "6264" 562157, \
 "9648" 562439, \
 "4694" 593171, \
 "7508" 593297, \
 "7138" 595024, \
 "6454" 595906, \
 "4676" 596233, \
 "5608" 596359, \
 "6650" 596551, \
 "5980" 596943, \
 "6146" 597175, \
 "8742" 597434, \
 "9530" 605281, \
 "5202" 627814, \
 "4054" 627971, \
 "4158" 628089, \
 "7926" 628210, \
 "7020" 631399, \
 "6644" 632088, \
 "7328" 632476, \
 "5784" 633688, \
 "6522" 633898, \
 "5006" 634246, \
 "9304" 634387, \
 "6318" 649744, \
 "7598" 650039, \
 "6392" 651969, \
 "6528" 652279, \
 "9586" 652628, \
 "4828" 677905, \
 "4040" 678036, \
 "3750" 678153, \
 "9662" 678262, \
 "6000" 710397, \
 "6930" 710631, \
 "7218" 711189, \
 "6464" 712195, \
 "4956" 712525, \
 "4316" 712664, \
 "9562" 712789, \
 "6182" 737214, \
 "5456" 737482, \
 "7446" 737664, \
 "6948" 739164, \
 "7936" 739740, \
 "8564" 742986, \
 "6654" 749439, \
 "5328" 749833, \
 "5402" 750001, \
 "6628" 750175, \
 "6942" 750557, \
 "5518" 751127, \
 "5354" 751313, \
 "5834" 751482, \
 "4484" 751698, \
 "4742" 751824, \
 "7310" 751950, \
 "5046" 753102, \
 "7126" 753246, \
 "4644" 754086, \
 "4382" 754212, \
 "4522" 754338, \
 "6822" 754464, \
 "6334" 754935, \
 "5220" 755233, \
 "7470" 755392, \
 "6828" 756973, \
 "7104" 757446, \
 "5972" 758256, \
 "9024" 758487, \
 "9282" 769370, \
 "8060" 784473, \
 "5898" 788185, \
 "8352" 788408, \
 "4112" 793525, \
 "4198" 793645, \
 "5616" 793767, \
 "6210" 793959, \
 "6374" 794232, \
 "6956" 794538, \
 "4228" 795134, \
 "6290" 795257, \
 "5104" 795547, \
 "4432" 795696, \
 "6484" 795822, \
 "5274" 796158, \
 "5520" 796321, \
 "5390" 796507, \
 "4366" 796680, \
 "" 797305 \
)
set size 1,1
set grid
set nokey
set border 10
set tics scale 0
set xlabel "Scer_chrI"
set ylabel "QRY"
set format "%.0f"
set xrange [1:230218]
set yrange [1:797305]
set linestyle 1  lt 1 lw 3 pt 6 ps 1
set linestyle 2  lt 3 lw 3 pt 6 ps 1
set linestyle 3  lt 2 lw 3 pt 6 ps 1
plot \
 "/Users/test/Desktop/UB/Advanced_Bioinformatics/exercise_02/blast/chrI-x-SRR6130428/chrI-x-SRR6130428_megablast.dnadiff.fplot" title "FWD" w lp ls 1, \
 "/Users/test/Desktop/UB/Advanced_Bioinformatics/exercise_02/blast/chrI-x-SRR6130428/chrI-x-SRR6130428_megablast.dnadiff.rplot" title "REV" w lp ls 2