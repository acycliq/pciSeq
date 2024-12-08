# pciSeq Diagnostics System

The pciSeq Diagnostics System implements a Model-View-Controller (MVC) architecture with Redis pub/sub for real-time communication between components. This design enables efficient monitoring and visualization of the pciSeq algorithm's diagnostic data.

## Architecture Overview

                          ┌─────────────────────┐
                          │    Redis Pub/Sub    │
                          │    Message Broker   │
                          └─────────┬───────────┘
                                    │
                         ┌──────────┼──────────┐
                         │                     │
                      Publish              Subscribe
                         ▲                     ▲
                         │                     │
                         │                     │
    ┌──────────┐    ┌─────────┐       ┌──────────────┐
    │          │    │         │       │              │
    │Algorithm ├───►│  Model  │       │    View      │
    │          │    │         │       │  Dashboard   │
    └──────────┘    └────┬────┘       └──────┬───────┘
                         ▲                   ▲
                         │    ┌──────────┐   │
                         └─── │Controller│───┘
                              │          │
                              └──────────┘



## Components

### Model (`model/diagnostic_model.py`)
- Takes algorithm data (gene efficiency and cell type data)
- Converts it to JSON strings with metadata
- Stores in Redis and notifies subscribers
- Initializes basic Redis connection

### View (`view/dashboard.py`)
- Implements Streamlit dashboard
- Subscribes to Redis channels
- Updates visualizations when new data arrives

### Controller (`controller/diagnostic_controller.py`)
Simple coordinator that:
- Launches the dashboard
- Passes algorithm updates to Model

### Redis
Separate database service that:
- Stores data as strings
- Enables pub/sub messaging between Model and View

## How It Works

1. **Startup**
   - Controller launches Streamlit dashboard
   - View subscribes to Redis channels
   - Redis server must be running

2. **During Algorithm Run**
   - pciSeq produces new data
   - Controller tells Model to publish it
   - Model converts to JSON and stores in Redis
   - Model notifies subscribers
   - View updates dashboard

## Redis Channels

Two main channels:
- `gene_efficiency`: Gene efficiency metrics
- `cell_type_posterior`: Cell type distribution data
