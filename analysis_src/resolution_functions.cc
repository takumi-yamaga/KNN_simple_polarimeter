// resolution_functions.cc

#include"resolution_functions.hh"

#include"TF1.h"
#include"TLorentzVector.h"
#include"TLatex.h"

#include"TRandom3.h"

// -----
Double_t cdc_resolution_phi_direction_rotated   = 1.75039e-1;
Double_t cdc_resolution_phirho_position_rotated = 1.01849e-2;
Double_t cdc_resolution_phi_rotation            = 1.42466e-1;
Double_t cdc_resolution_z_direction_rotated     = 7.98753e-1;
Double_t cdc_resolution_z_position_rotated      = 1.06003e-1;
Double_t cdc_resolution_z_rotation              = 4.06320e-1;
// -----
Double_t tracker_resolution_phi_position = 1.; // mm (sigma)
Double_t tracker_resolution_z_position   = 1.; // mm (sigma)
// -----

void CDCResolution(TVector3& pos, TVector3& dir){
    if(pos.Mag()==0||dir.Mag()==0){
        return;
    }

    TVector3 pos_org = pos; 
    TVector3 dir_org = dir; 

    TRandom3* random = new TRandom3(0);

    //// phi ---------------------
    double diff_phi_direction_rotated = random->Gaus(0.,cdc_resolution_phi_direction_rotated);
    double diff_phirho_position_rotated = random->Gaus(0.,cdc_resolution_phirho_position_rotated);
    double rotation_phi_vs_phirho = atan(cdc_resolution_phi_rotation);

    double diff_phi_direction = cos(rotation_phi_vs_phirho) * diff_phi_direction_rotated - sin(rotation_phi_vs_phirho) * diff_phirho_position_rotated;
    double diff_phirho_position = sin(rotation_phi_vs_phirho) * diff_phi_direction_rotated + cos(rotation_phi_vs_phirho) * diff_phirho_position_rotated;

    double phi_direction = dir_org.Phi() + diff_phi_direction*TMath::Pi()/180.;
    double phirho_position = pos_org.Phi() * pos_org.Perp() + diff_phirho_position;
    //// -------------------------

    //// z -----------------------
    double diff_theta_direction_rotated = random->Gaus(0.,cdc_resolution_z_direction_rotated);
    double diff_z_position_rotated = random->Gaus(0.,cdc_resolution_z_position_rotated);
    double rotation_theta_vs_z = atan(-cdc_resolution_z_rotation);

    double diff_theta_direction = cos(rotation_theta_vs_z) * diff_theta_direction_rotated - sin(rotation_theta_vs_z) * diff_z_position_rotated;
    double diff_z_position = sin(rotation_theta_vs_z) * diff_theta_direction_rotated + cos(rotation_theta_vs_z) * diff_z_position_rotated;

    double theta_direction = dir_org.Theta() + diff_theta_direction*TMath::Pi()/180.;
    double z_position = pos_org.Z() + diff_z_position;
    //// -------------------------

    double dir_x = dir_org.Mag() * sin(theta_direction) * cos(phi_direction);
    double dir_y = dir_org.Mag() * sin(theta_direction) * sin(phi_direction);
    double dir_z = dir_org.Mag() * cos(theta_direction);
    dir.SetXYZ(dir_x,dir_y,dir_z);

    double pos_x = pos_org.Perp() * cos(phirho_position/pos_org.Perp());
    double pos_y = pos_org.Perp() * sin(phirho_position/pos_org.Perp());
    double pos_z = z_position;
    pos.SetXYZ(pos_x,pos_y,pos_z);

    delete random;

    return;
}

void TrackerResolution(TVector3& pos) {
    if(pos.Mag()==0){
        return;
    }
    TVector3 pos_org = pos;

    TRandom3* random = new TRandom3(0);
    double diff_phirho_dir = random->Gaus(0.,tracker_resolution_phi_position);
    double phirho_position = pos_org.Phi()*pos_org.Perp() + diff_phirho_dir;
    double diff_z_position = random->Gaus(0.,tracker_resolution_z_position);
    double z_position = pos_org.Z() + diff_z_position;

    double pos_x = pos_org.Perp() * cos(phirho_position/pos_org.Perp());
    double pos_y = pos_org.Perp() * sin(phirho_position/pos_org.Perp());
    double pos_z = z_position;
    pos.SetXYZ(pos_x,pos_y,pos_z);

    delete random;

    return;
}
