# pciSeq Diagnostics System

This folder contains the implementation of the diagnostics system for the pciSeq project. The system follows the Model-View-Controller (MVC) pattern and uses Redis for real-time data communication.

## Structure

1. Model: `model.py`
2. View: `diagnostics.py`
3. Controller: `launch_diagnostics.py`
4. Data Publisher: `redis_publisher.py`
5. Utilities: `utils.py`
6. Configuration: `config.py`

## Components

### Model (`model.py`)
- Manages data and business logic related to diagnostics
- Interacts with Redis database to store and retrieve diagnostic information
- Provides methods for accessing data

### View (`diagnostics.py`)
- Implements a Streamlit dashboard for data visualization
- Subscribes to Redis channels for real-time updates
- Displays gene efficiency and cell type distribution charts

### Controller (`launch_diagnostics.py`)
- Coordinates between the Model and View
- Handles the launching and shutting down of the dashboard

### Data Publisher (`redis_publisher.py`)
- Publishes diagnostic data to Redis
- Acts as part of the Controller in the MVC pattern
- Updates the Model (Redis database) with the latest algorithm state

### Utilities (`utils.py`)
- Provides helper functions for Redis operations
- Handles Redis installation, starting, and stopping across different platforms

### Configuration (`config.py`)
- Contains settings for Redis connection

## Data Flow
1. `redis_publisher.py` publishes diagnostic data to Redis
2. `model.py` retrieves data from Redis and provides it to the View
3. `diagnostics.py` (View) subscribes to Redis channels and updates the dashboard in real-time

## Usage
The diagnostics system is typically launched as part of the main pciSeq algorithm. It can also be run independently for development or testing purposes.

To launch the dashboard manually:

`from pciSeq.src.diagnostics.launch_diagnostics import launch_dashboard`

`launch_dashboard()`

Note: Ensure that Redis is installed and running before launching the diagnostics dashboard.

## Redis Requirements
The diagnostics system requires Redis to be installed and running. The `utils.py` file provides functions to check Redis status, install Redis if necessary, and start/stop the Redis server.

For Windows users, the system uses Memurai, a Redis-compatible server. The installation file is included in the project.

For Linux and macOS users, the system attempts to use the native Redis server.

## Troubleshooting
If you encounter issues with the diagnostics dashboard:
1. Check if Redis is installed and running using the functions in `utils.py`
2. Ensure that the Redis connection settings in `config.py` are correct
3. Check the logs for any error messages related to Redis or Streamlit

For more detailed information about each component, refer to the docstrings and comments within the individual files.


