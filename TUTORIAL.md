# CoLoMoTo Docker Tutorial (macOS)

## Overview

The CoLoMoTo Docker provides a Jupyter environment with pre-installed tools for boolean network analysis (ginsim, biolqm, boolsim, etc.). This guide covers usage on macOS.

## Prerequisites

Make sure Docker Desktop for Mac is installed and running:

- Download from: <https://www.docker.com/products/docker-desktop>
- After installation, ensure Docker Desktop is running (check the menu bar for the Docker whale icon)

## Installation Methods

### Method 1: Python Helper Script (Recommended)

Install the colomoto-docker helper tool:

```bash
pip install -U colomoto-docker
```

**To start from the RookFields root directory (recommended):**

```bash
colomoto-docker -v notebooks:local-notebooks
```

This command:
- Mounts your `notebooks/` directory as `local-notebooks` in the container
- Opens JupyterLab in your default browser
- Uses the latest CoLoMoTo Docker image

Your notebooks will appear in the `local-notebooks` directory in JupyterLab.

**Useful Options:**
- `-v notebooks:local-notebooks` - Mount local notebooks directory
- `--lab` - Use JupyterLab interface (default)
- `--notebook` - Use classic Jupyter Notebook interface
- `-V 2025-03-01` - Specify a specific image version
- `--shell` - Open a command shell instead

### Method 2: Direct Docker Commands

If you prefer not using the Python helper, use Docker directly:

```bash
# Pull the image first (one time)
docker pull colomoto/colomoto-docker:2025-03-01

# Run from the RookFields root directory
docker run -it --rm -p 8888:8888 -v "$(pwd)/notebooks":/notebook colomoto/colomoto-docker:2025-03-01
```

**Docker command breakdown:**

- `-it`: Interactive terminal
- `--rm`: Remove container when stopped (cleanup)
- `-p 8888:8888`: Map port 8888 (Jupyter default)
- `-v "$(pwd)/notebooks":/notebook`: Mount your local `notebooks/` directory to `/notebook` in the container
- `colomoto/colomoto-docker:2025-03-01`: The CoLoMoTo Docker image with specific version tag

**Important notes:**

- **Quotes around paths:** The quotes in `"$(pwd)/notebooks"` handle spaces in macOS paths
- **Image versions:** Check [Docker Hub](https://hub.docker.com/r/colomoto/colomoto-docker/tags) for available tags
- **macOS compatibility:** On older macOS systems, you may need to add `--ulimit nofile=8096` if you encounter file descriptor issues

## Accessing Jupyter

After starting the container, navigate to:

```
http://localhost:8888
```

or

```
http://127.0.0.1:8888/lab
```

in your web browser (Safari, Chrome, etc.).

**Important:** The CoLoMoTo Docker runs without authentication by default. No token is required.

**macOS Note:** If you're using older Docker Toolbox instead of Docker Desktop, use `docker-machine ip default` to find the correct IP address instead of localhost.

## Opening Your Notebooks

Once in JupyterLab:

1. Navigate to the file browser on the left
2. Your notebooks from `notebooks/` will be available (in `/notebook` if using direct Docker, or in the current directory if using the Python helper)
3. Open `boolean_dsgrn-Copy1.ipynb` or any other notebook

**Kernel Issues:**

If the kernel fails to start or shows continuous "Nudge" warnings and TimeoutError:

1. **Try creating a new notebook first:**
   - In JupyterLab, create a new notebook (File → New → Notebook)
   - Select the Python 3 kernel
   - Try running a simple cell: `print("test")`
   - If this works, the issue is with your existing notebook

2. **If new notebooks also fail:**
   - Stop the container (Ctrl+C)
   - Clear any stale containers: `docker ps -a` and `docker rm <container-id>` if needed
   - Try using the classic notebook interface instead: `colomoto-docker --notebook -v notebooks:local-notebooks`
   - Try a different image version: `colomoto-docker -V 2024-05-01 -v notebooks:local-notebooks`

3. **If the specific notebook is the problem:**
   - The notebook may have corrupted kernel metadata
   - Try opening it in a text editor and check for Python syntax errors in code cells
   - Create a fresh notebook and copy cells over one by one

4. **Docker resource limits:**
   - Check Docker Desktop → Settings → Resources
   - Increase memory allocation if low (recommend at least 4GB)
   - Restart Docker Desktop after changing settings

## Stopping CoLoMoTo Docker

- **With Python helper:** Press `Ctrl+C` in the terminal
- **With direct Docker:** Press `Ctrl+C` twice to stop the server

The `--rm` flag in Docker commands ensures the container is automatically removed on exit.

## Available Tools in CoLoMoTo

- **ginsim**: Gene regulatory network visualization and analysis
- **biolqm**: Logical qualitative modeling toolkit
- **boolsim**: Boolean network simulator
- **pint**: Process hitting temporal analysis
- And many more CoLoMoTo tools

## Important Notes

### File Persistence

⚠️ **Critical:** Files within the container are deleted after stopping, **except** those in mounted directories.

- **Persists:** Files in your mounted `notebooks/` directory (or wherever you mounted with `-v`)
- **Does not persist:** Files created elsewhere in the container (unless you use the `persistent` directory inside the container)

Always save your work in the mounted directory to avoid losing it.

### Python Environment

- All Python packages depending on CoLoMoTo tools (ginsim, biolqm, boolsim) only work inside the CoLoMoTo Docker environment
- The Docker container has its own Python environment, separate from your local Python installation
- Import statements like `import ginsim` will fail outside the container

### macOS Compatibility

- **Apple Silicon (M1/M2/M3):** Docker runs the container using Rosetta 2 translation (x86_64 emulation) - this is normal and works fine
- **Firewall:** You may need to allow Docker network access through the firewall if prompted
- **Performance:** First-time kernel startup may be slow; subsequent starts are faster
