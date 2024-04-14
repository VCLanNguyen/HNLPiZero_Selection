col_dict = {
        'Teal':("#005A8F")
        , 'Aqua':("#00BBCC")
        , 'SkyBlue':("#64CBE8")
        , 'MintGreen':("#00CE7C")
        , 'Spearmint':("#3BD4AE")
        , 'PastelGreen':("#A1DED2")
        , 'Mauve':("#663DB3")
        , 'Purple':("#981F92")
        , 'Lavender':("#DAA8E2")
        , 'Coral':("#E7004C")
        , 'Flamingo':("#FF6371")
        , 'Peach':("#FF9664")
        , 'DeepViolet':("#440099")
        , 'PowderBlue':("#9ADBE8")
        , 'MidnightBlack':("#131E29")
        , 'RosyBrown4':("#8B6969")
        , 'SlateGray':("#708090")
        , 'LightGray':("#D0D2D4")
    }

event_col_dict = {                
                0:  col_dict['Flamingo']       
                , 1: col_dict['Coral']         
                , 2: col_dict['Mauve']         
                , 3: col_dict['Teal']        
                , 4: col_dict['PowderBlue']        
                , 5: col_dict['MintGreen']        
                , 6: col_dict['Lavender']        
                , 7: col_dict['Peach']        
                , 8: col_dict['RosyBrown4']        
                , 9: col_dict['LightGray']        
                , -1: col_dict['MidnightBlack'] 
                , -99: col_dict['DeepViolet'] 
            }

event_label_dict = {                
                0:  "HNL ${\pi}^{0}$"       
                , 1:  "Non-FV HNL ${\pi}^{0}$"        
                , 2:  "Dirt HNL"         
                , 3:  "NC ${\pi}^{0}$"        
                , 4:  "Other NC"         
                , 5:  "CC $\\nu_{\mu}$"          
                , 6:  "CC $\\nu_{e}$"           
                , 7:  "Non-FV $\\nu$"         
                , 8:  "Dirt $\\nu$"          
                , 9:  "Cosmic"          
                , -1: "Unknown"
                , -99: "Bad Reco Signal"
            }

event_type = [ 0
              #, 1
              #, 2
              , 3
              , 4
              , 5
              , 6
              , 7
              , 8
              , 9
              #, -99
              #, -1
              ]

event_col = [
                event_col_dict[event_type[0]]
                , event_col_dict[event_type[1]]
                , event_col_dict[event_type[2]]
                , event_col_dict[event_type[3]]
                , event_col_dict[event_type[4]]
                , event_col_dict[event_type[5]]
                , event_col_dict[event_type[6]]
                , event_col_dict[event_type[7]]
                #, event_col_dict[event_type[8]]
                #, event_col_dict[event_type[9]]
                #, event_col_dict[event_type[10]]
            ]
            
event_label = [
                event_label_dict[event_type[0]]
                , event_label_dict[event_type[1]]
                , event_label_dict[event_type[2]]
                , event_label_dict[event_type[3]]
                , event_label_dict[event_type[4]]
                , event_label_dict[event_type[5]]
                , event_label_dict[event_type[6]]
                , event_label_dict[event_type[7]]
                #, event_label_dict[event_type[8]]
                #, event_label_dict[event_type[9]]
                #, event_label_dict[event_type[10]]
            ]





