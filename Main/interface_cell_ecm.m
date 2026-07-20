function C_interface = interface_cell_ecm(C_cell_end, C_ecm_start, D_cell, D_ecm)
%% Fick's-law flux-continuity value at the cell-ECM interface.
%
%   Inputs:
%     C_cell_end  : last interior point of the cell region  (adjacent to interface)
%     C_ecm_start : first interior point of the ECM region  (adjacent to interface)
%     D_cell      : diffusion coefficient in the cell region
%     D_ecm       : diffusion coefficient in the ECM region
%
%   Output:
%     C_interface : concentration to impose as the shared interface boundary
%                   (BC_R for the cell, BC_L for the ECM)
%
%%
    C_interface = (D_cell*C_cell_end + D_ecm*C_ecm_start) / (D_cell + D_ecm);

end