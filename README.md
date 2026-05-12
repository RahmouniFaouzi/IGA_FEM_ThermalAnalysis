#  Heat Transfer in Composite Materials Using FEM, Isogeometric Analysis, and Symbolic Regression

This repository contains the MATLAB code developed for the publication:

**Comparative Study of Steady-State and Transient Heat Transfer in Composites Using FEM, Isogeometric Analysis, and Physics-Informed Machine Learning Correlation**

**Authors:** Dr. Faouzi Rahmouni, Prof. Amar Khennane

<img width="280" height="204" alt="image" src="https://github.com/user-attachments/assets/435dd904-8b29-4fc7-9737-0e6b2af5dbb0" />
<img width="280" height="204" alt="image" src="https://github.com/user-attachments/assets/d0d6f266-2ea4-45fe-85e8-b815cb61952f" />
<img width="280" height="204" alt="image" src="https://github.com/user-attachments/assets/e37be941-62e7-4a67-bf37-8d7685dfec66" />
<img width="280" height="204" alt="image" src="https://github.com/user-attachments/assets/5debf457-9241-4567-8018-3b41514601f9" />
<img width="240" height="190" alt="image" src="https://github.com/user-attachments/assets/53d24d8f-059d-40f0-88dd-64f030d79e15" />
<img width="240" height="204" alt="image" src="https://github.com/user-attachments/assets/22ad4c49-16ca-428b-93c7-e1ae897b4d35" />
<img width="674" height="363" alt="image" src="https://github.com/user-attachments/assets/44c2fd9f-4054-48f4-8697-da8db3f43e56" />

## Summary
This code supports the publication and reproduces the main numerical results for steady-state and transient heat transfer in composite plates using:

- **FEM**
- **Isogeometric Analysis (IGA)**:
  - NURBS
  - HB-Splines
  - Graded-NURBS
- **Physics-Informed Machine Learning Correlation**

## Main Applications
- Isotropic, orthotropic, and laminated MMC plates
- Plates with and without a central hole
- Prescribed temperature, insulation, and mixed conduction-convection boundary conditions

## How to Use

Run the MATLAB scripts in the relevant folder depending on the method and case study.

### Examples
- `FEM/` for finite element benchmarks
- `IGA/` for NURBS, HB-Spline, and holed-plate models
- `Machine Learning/` for Physics-Informed Machine Learning Correlation

## Purpose

This repository is provided to:
- Support reproducibility of the publication
- Compare FEM and multi-formulation IGA
- Provide benchmark MATLAB examples for thermal analysis in composite materials

## Citation

If you use this code, please cite the associated publication.

## License

This code is provided for academic and research purposes.
