
        from pymol import cmd,stored
        
        set depth_cue, 1
        set fog_start, 0.4
        
        set_color b_col, [36,36,85]
        set_color t_col, [10,10,10]
        set bg_rgb_bottom, b_col
        set bg_rgb_top, t_col      
        set bg_gradient
        
        set  spec_power  =  200
        set  spec_refl   =  0
        
        load "data/8WCC_fixed.pdb", protein
        create ligands, protein and organic
        select xlig, protein and organic
        delete xlig
        
        hide everything, all
        
        color white, elem c
        color bluewhite, protein
        #show_as cartoon, protein
        show surface, protein
        #set transparency, 0.15
        
        show sticks, ligands
        set stick_color, magenta
        
        
        
        
        # SAS points
 
        load "data/8WCC_fixed.pdb_points.pdb.gz", points
        hide nonbonded, points
        show nb_spheres, points
        set sphere_scale, 0.2, points
        cmd.spectrum("b", "green_red", selection="points", minimum=0, maximum=0.7)
        
        
        stored.list=[]
        cmd.iterate("(resn STP)","stored.list.append(resi)")    # read info about residues STP
        lastSTP=stored.list[-1] # get the index of the last residue
        hide lines, resn STP
        
        cmd.select("rest", "resn STP and resi 0")
        
        for my_index in range(1,int(lastSTP)+1): cmd.select("pocket"+str(my_index), "resn STP and resi "+str(my_index))
        for my_index in range(1,int(lastSTP)+1): cmd.show("spheres","pocket"+str(my_index))
        for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_scale","0.4","pocket"+str(my_index))
        for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_transparency","0.1","pocket"+str(my_index))
        
        
        
        set_color pcol1 = [0.361,0.576,0.902]
select surf_pocket1, protein and id [4165,4228,4232,1650,1653,4206,4210,4212,2963,3073,4640,4644,1590,1636,1639,1640,1641,1258,1261,4556,4560,1588,2919,1657,1593,1597,1601,4564,4274,1313,4559,1708,1706,1703,1647,1644,1710,1755,3134,3137,2466] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.349,0.702]
select surf_pocket2, protein and id [1478,1481,2825,2780,2822,1410,1534,2880,2881,2882,2887,2889,2762,2888,2877,1405,1486,2623,2628,2754,2758,2905,2640,2646,2592,2597,2637,2748,2749,2912,1477,1545,1548,1551,1597,1599,1601,2536,2530] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.404,0.361,0.902]
select surf_pocket3, protein and id [4319,4395,4399,4390,4391,4455,4456,4427,4433,4444,4445,4324,2955,2959,2906,2909,2915,4208,4272,4509,4570,4572,2806,2813,2812,2814,2935,2927,2945,2950,2951,4269,4333,4266,4471,1370,1374,2871,4459,4460,2894,1325,4519,1305,1308,4524,1367,1301,1355,2862,2805,2884,2888] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.416,0.278,0.702]
select surf_pocket4, protein and id [1178,1183,1187,1622,1624,1626,1206,1566,1177,1454,1455,1457,1515,1521,1512,1418,1511,1580,1583,1574,1576,1459,1461,1463,1247,1250,1586] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.663,0.361,0.902]
select surf_pocket5, protein and id [4905,4907,910,912,914,1905,2086,2084,1968,1921,1923,4870,4786,4787,4793,4864,4868,986,990,1853,1856] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.616,0.278,0.702]
select surf_pocket6, protein and id [1918,3337,4797,4799,4795,4029,4028,4730,4085,3335,3259,1856] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.361,0.878]
select surf_pocket7, protein and id [1665,1666,1670,1687,2430,2424,2499,2497,1183,1624,1626,1689,1693,2361,2372,2376,1125,1127] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.278,0.584]
select surf_pocket8, protein and id [1915,1918,1921,1923,4870,3987,3993,3400,3396,4866,3974,3970,3924,3470,3928,3989,1984] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.902,0.361,0.620]
select surf_pocket9, protein and id [509,4591,443,456,4584,4585,4587,4658,4580,4653,4581,4530,4609,4613,4531,4527,460,386,388,392,1328,4519,4522,4523] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.278,0.380]
select surf_pocket10, protein and id [3323,3324,3340,3342,3346,3350,3354,3984,3325,3410,3413,3414,3416,3418,3422,3420,3439,3443,3989] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.902,0.361,0.361]
select surf_pocket11, protein and id [1890,1839,1843,2337,2402,3239,1821,1768,3200,3207,1771,1775,1777,3237,3243,3185,3191] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.702,0.380,0.278]
select surf_pocket12, protein and id [2395,2397,2396,2456,2459,2460,2462,2464,1775,2410,2421,1717,1777,2411,3120,3137,3122,3124,2402,2406] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.902,0.620,0.361]
select surf_pocket13, protein and id [2231,2254,2306,944,948,937,938,941,945,947,950,2191,2184,952,1004,1008,762] 
set surface_color,  pcol13, surf_pocket13 
set_color pcol14 = [0.702,0.584,0.278]
select surf_pocket14, protein and id [3155,3157,4114,3151,3169,3173,3091,3093,4249,3219,3218,3229,3215,3144,3275,3279] 
set surface_color,  pcol14, surf_pocket14 
set_color pcol15 = [0.902,0.878,0.361]
select surf_pocket15, protein and id [4900,749,4923,4808,986] 
set surface_color,  pcol15, surf_pocket15 
   
        
        deselect
        
        orient
        