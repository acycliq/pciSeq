# pciSeq Diagnostics System

This folder contains the implementation of the diagnostics system for the pciSeq project, following the Model-View-Controller (MVC) pattern with Redis for real-time data communication.

## Components

### Model (`model/diagnostic_model.py`)
- Manages Redis connection and data operations
- Publishes diagnostic data to Redis channels
- Retrieves data for the View component

### View (`view/dashboard.py`)
- Implements Streamlit dashboard
- Subscribes to Redis for real-time updates
- Creates interactive visualizations using Altair

### Controller (`controller/diagnostic_controller.py`)
- Coordinates Model and View components
- Manages dashboard lifecycle
- Handles system signals for clean shutdown
- Updates diagnostic data from algorithm


