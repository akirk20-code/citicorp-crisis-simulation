# Citicorp Center CFD Simulation

A Computational Fluid Dynamics (CFD) project simulating atmospheric boundary layer flow around the Citicorp Center in New York City. This project uses OpenFOAM (ESI version) for the simulation and Python for procedural geometry generation using NYC Open Data.

## Prerequisites

*   **OpenFOAM:** v2406 (ESI) or compatible (requires `foamRun`, `incompressibleFluid` solver).
*   **Python:** 3.8+
*   **System:** Linux or WSL2 (Windows Subsystem for Linux).

## Quick Start

1.  **Install Python Dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

2.  **Run Simulation:**
    The `Allrun` script handles geometry generation, meshing, and solving.
    ```bash
    ./Allrun
    ```

3.  **Clean Case:**
    To remove all generated mesh and solution files:
    ```bash
    ./Allclean
    ```

## Project Structure

*   `generate_stl.py`: Fetches building footprints from NYC Open Data and generates `.stl` files.
*   `system/`: OpenFOAM control dictionaries (`controlDict`, `snappyHexMeshDict`, etc.).
*   `constant/`: Physical properties and turbulence models.
*   `0/`: Initial boundary conditions.

## Simulation Details

*   **Solver:** `foamRun` (incompressible, steady-state).
*   **Turbulence Model:** k-omega SST.
*   **Mesh:** `snappyHexMesh` with background `blockMesh`.
*   **Geometry:** Citicorp tower (detailed) + surrounding buildings (extruded footprints).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
