#events tree
mct_branches = [
                "run"
                , "subrun"
                , "event"
                ,"nu_mctruth_id"
                ,"nu_event_type"
                ,"nu_en_dep"
            ]

#mevprtl
mevprtl_branches =[
                    "run"
                    , "subrun"
                    , "event"
                    , "mevprtl_decay_pos_x"
                    , "mevprtl_decay_pos_y"
                    , "mevprtl_decay_pos_z"
                    , "mevprtl_decay_pos_t"
                    , "mevprtl_flux_weight"
                    , "mevprtl_ray_weight"
                    , "mevprtl_decay_weight"
                ]

#subruns tree
subrun_branches = [
                    "run"
                    , "subrun"
                    , "pot"
                    , "spills"
                    , "ngenevts"
                ]

#fluxsys
fluxw_branches = [
                "run",
                "subrun",
                "event",
                "slc_flux_weight_expskin"
                "slc_flux_weight_horncurrent"
                "slc_flux_weight_kminus"
                "slc_flux_weight_kplus"
                "slc_flux_weight_kzero"
                "slc_flux_weight_nucleoninexsec"
                "slc_flux_weight_nucleonqexsec"
                "slc_flux_weight_nucleontotxsec"
                "slc_flux_weight_piminus"
                "slc_flux_weight_pioninexsex"
                "slc_flux_weight_pionqexsec"
                "slc_flux_weight_piontotxsec"
                "slc_flux_weight_piplus"
                "slc_flux_weight_total"
                "slc_xsec_unisim_DecayAngMEC"
                "slc_xsec_unisim_ThetaDelta2NRad"
                "slc_xsec_unisim_Theta_Delta2Npi"
                "slc_xsec_unisim_VecFFCCQEshape"
                "slc_xsec_multisigma_CoulombCCQE"
                "slc_xsec_multisigma_NonRESBGvbarnCC1pi"
                "slc_xsec_multisigma_NonRESBGvbarnCC2pi"
                "slc_xsec_multisigma_NonRESBGvbarnNC1pi"
                "slc_xsec_multisigma_NonRESBGvbarnNC2pi"
                "slc_xsec_multisigma_NonRESBGvbarpCC1pi"
                "slc_xsec_multisigma_NonRESBGvbarpCC2pi"
                "slc_xsec_multisigma_NonRESBGvbarpNC1pi"
                "slc_xsec_multisigma_NonRESBGvbarpNC2pi"
                "slc_xsec_multisigma_NonRESBGvnCC1pi"
                "slc_xsec_multisigma_NonRESBGvnCC2pi"
                "slc_xsec_multisigma_NonRESBGvnNC1pi"
                "slc_xsec_multisigma_NonRESBGvnNC2pi"
                "slc_xsec_multisigma_NonRESBGvpCC1pi"
                "slc_xsec_multisigma_NonRESBGvpCC2pi"
                "slc_xsec_multisigma_NonRESBGvpNC1pi"
                "slc_xsec_multisigma_NonRESBGvpNC2pi"
                "slc_xsec_multisigma_NormCCMEC"
                "slc_xsec_multisigma_NormNCMEC"
                "slc_xsec_multisigma_RDecBR1eta"
                "slc_xsec_multisigma_RDecBR1gamma"
                "slc_xsec_multisigma_RPA_CCQE"
                "slc_xsec_multisigma_NormNCCOH"
                "slc_xsec_multisigma_NormCCCOH"
                "slc_xsec_multisim_ZExpA"
                "slc_xsec_multisim_NCEL"
                "slc_xsec_multisim_CCRES"
                "slc_xsec_multisim_NCRES"
                "slc_xsec_multisim_DISBY"
                "slc_xsec_multisim_FSI_pi"
                "slc_xsec_multisim_FSI_N"
                "slc_xsec_multisim_total"
                "slc_geant4_multisim_reinteractions"
                ]

#slice
slc_branches = [
                "run"
                , "subrun"
                , "event"
                , "slc_id"
                #, "n_slc"
                , "slc_n_pfps"
                , "slc_is_clear_cosmics"                                         
                #, "slc_primary_pfp_id"                                           
                #, "slc_primary_pfp_pdg"                                          
                #, "slc_n_primary_daughters"                                      
                , "slc_vtx_x"                                                    
                , "slc_vtx_y"                                                    
                , "slc_vtx_z"                                                    
                , "slc_is_fv"                                                    
                , "slc_crumbs_score"                                             
                #, "slc_crumbs_nc_score"                                          
                #, "slc_crumbs_ccnue_score"                                       
                #, "slc_opt0_time"                                                
                , "slc_opt0_score"                                               
                , "slc_opt0_measPE"                                              
                #, "slc_opt0_hypoPE"         
                , "slc_opt0_frac"         
                , "slc_opt0_time_corrected_Z_pandora"
                #, "slc_n_trks"                                                   
                #, "slc_n_shws"                 
                #, "slc_n_stub"
                , "slc_total_shower_E"
                , "slc_total_track_E"
                #, "slc_n_primary_trks"                                           
                #, "slc_n_primary_shws"                                           
                #, "slc_n_dazzle_muons"                                           
                #, "slc_n_dazzle_pions"                                           
                #, "slc_n_dazzle_pions_thresh"                                    
                #, "slc_n_dazzle_protons"                                         
                #, "slc_n_dazzle_protons_thresh"                                  
                #, "slc_n_dazzle_other"                                           
                #, "slc_n_primary_dazzle_muons"                                   
                #, "slc_n_primary_dazzle_pions"                                   
                #, "slc_n_primary_dazzle_pions_thresh"                            
                #, "slc_n_primary_dazzle_protons"                                 
                #, "slc_n_primary_dazzle_protons_thresh"                          
                #, "slc_n_primary_dazzle_other"                                   
                #, "slc_n_razzle_electrons"                                       
                #, "slc_n_razzle_photons"                                         
                #, "slc_n_razzle_other"                                           
                #, "slc_n_primary_razzle_electrons"                               
                #, "slc_n_primary_razzle_photons"                                 
                #, "slc_n_primary_razzle_other"                                   
                , "slc_n_razzled_electrons"                                      
                , "slc_n_razzled_muons"                                          
                , "slc_n_razzled_photons"                                        
                #, "slc_n_razzled_pions"                                          
                , "slc_n_razzled_pions_thresh"                                   
                #, "slc_n_razzled_protons"                                        
                , "slc_n_razzled_protons_thresh"
                #, "slc_n_primary_razzled_electrons"     
                #, "slc_n_primary_razzled_muons"         
                #, "slc_n_primary_razzled_photons"       
                #, "slc_n_primary_razzled_pions"         
                #, "slc_n_primary_razzled_pions_thresh"  
                #, "slc_n_primary_razzled_protons"       
                #, "slc_n_primary_razzled_protons_thresh"
                #TRUE
                , "slc_comp"           
                #, "slc_pur"            
                #, "slc_true_mctruth_id"
                , "slc_true_event_type"
                #, "slc_true_en_dep"    
                #, "slc_true_vtx_x"     
                #, "slc_true_vtx_y"     
                #, "slc_true_vtx_z"     
                #, "slc_true_vtx_t"     
                #, "slc_true_vtx_t"     
                #, "slc_true_vtx_t_corrected_Z"
                ]

#pfp
pfp_branches = [
                "run"
                , "subrun"
                , "event"
                #, "slc_pfp_id"                              
                #, "slc_pfp_pdg"                             
                #, "slc_pfp_n_children"                      
                #, "slc_pfp_primary"                         
                #, "slc_pfp_primary_child"                   
                , "slc_pfp_n_hits"                          
                #, "slc_pfp_n_sps"                           
                , "slc_pfp_track_score"                     
                #, "slc_pfp_good_track"                      
                #, "slc_pfp_good_shower"                     
                #, "slc_pfp_comp"                            
                #, "slc_pfp_pur"                             
                #, "slc_pfp_cnnscore_track"                  
                #, "slc_pfp_cnnscore_shower"                 
                #, "slc_pfp_cnnscore_noise"                  
                #, "slc_pfp_cnnscore_michel"                 
                #, "slc_pfp_cnnscore_endmichel"              
                #, "slc_pfp_cnnscore_nclusters"              
                , "slc_pfp_razzled_electron_score"          
                , "slc_pfp_razzled_muon_score"              
                , "slc_pfp_razzled_photon_score"            
                , "slc_pfp_razzled_pion_score"              
                , "slc_pfp_razzled_proton_score"            
                #, "slc_pfp_razzled_pdg"                     
                #TRUE                                                
                #, "slc_pfp_true_trackid"                    
                #, "slc_pfp_true_pdg"                        
                #, "slc_pfp_true_energy"                     
                #, "slc_pfp_true_p_x"                        
                #, "slc_pfp_true_p_y"                        
                #, "slc_pfp_true_p_z"                        
                #TRACK
                #, "slc_pfp_track_start_x"                  
                #, "slc_pfp_track_start_y"                  
                #, "slc_pfp_track_start_z"                  
                #, "slc_pfp_track_end_x"                    
                #, "slc_pfp_track_end_y"                    
                #, "slc_pfp_track_end_z"                    
                #, "slc_pfp_track_dir_x"                    
                #, "slc_pfp_track_dir_y"                    
                #, "slc_pfp_track_dir_z"                    
                , "slc_pfp_track_length"                   
                #, "slc_pfp_track_dazzle_muon_score"        
                #, "slc_pfp_track_dazzle_pion_score"        
                #, "slc_pfp_track_dazzle_proton_score"      
                #, "slc_pfp_track_dazzle_other_score"       
                #, "slc_pfp_track_dazzle_pdg"               
                , "slc_pfp_track_ke"                       
                #, "slc_pfp_track_charge"                   
                #, "slc_pfp_track_chi2_muon"                
                #, "slc_pfp_track_chi2_pion"                
                #, "slc_pfp_track_chi2_kaon"                
                #, "slc_pfp_track_chi2_muon"                
                #, "slc_pfp_track_chi2_pdg"                 
                #, "slc_pfp_track_mcs_mean_scatter"         
                #, "slc_pfp_track_mcs_max_scatter_ratio"    
                #, "slc_pfp_track_range_p"                  
                #, "slc_pfp_track_closest_approach_mean_dca"
                #, "slc_pfp_track_stopping_dedx_chi2_ratio" 
                #, "slc_pfp_track_stopping_dedx_pol0_fit"   
                , "slc_pfp_track_theta"
                , "slc_pfp_track_phi"
                , "slc_pfp_track_contained"
                #SHOWER
                #, "slc_pfp_shower_start_x"              
                #, "slc_pfp_shower_start_y"              
                #, "slc_pfp_shower_start_z"              
                #, "slc_pfp_shower_end_x"                
                #, "slc_pfp_shower_end_y"                
                #, "slc_pfp_shower_end_z"                
                , "slc_pfp_shower_conv_gap"             
                #, "slc_pfp_shower_dir_x"                
                #, "slc_pfp_shower_dir_y"                
                #, "slc_pfp_shower_dir_z"                
                , "slc_pfp_shower_length"               
                , "slc_pfp_shower_open_angle"           
                , "slc_pfp_shower_energy"               
                , "slc_pfp_shower_dedx"   
                #, "slc_pfp_shower_sqrt_energy_density"  
                #, "slc_pfp_shower_modified_hit_density" 
                #, "slc_pfp_shower_razzle_electron_score"
                #, "slc_pfp_shower_razzle_photon_score"  
                #, "slc_pfp_shower_razzle_other_score"   
                #, "slc_pfp_shower_razzle_pdg"           
                #, "slc_pfp_shower_cosmic_dist"          
                #, "slc_pfp_shower_track_length"         
                #, "slc_pfp_shower_track_width"          
                #, "slc_pfp_shower_density_grad"         
                #, "slc_pfp_shower_density_pow"          
                , "slc_pfp_shower_theta"
                , "slc_pfp_shower_phi"
                , "slc_pfp_shower_contained"
                ]

