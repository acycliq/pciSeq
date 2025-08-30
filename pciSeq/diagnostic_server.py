"""
pciSeq Diagnostic Server

Provides HTTP API endpoints for running check_cell and check_spot diagnostics.

Usage:
  import pciSeq
  pciSeq.diagnostic_server('/path/to/varbayes.joblib')
"""

from flask import Flask, request, jsonify
from flask_cors import CORS
import joblib
import pickle
import os
import sys
import base64
import io
import logging

# Configure matplotlib for headless operation
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Set up logging
logger = logging.getLogger(__name__)

app = Flask(__name__)
CORS(app)  # Enable cross-origin requests

# Global variable for the VarBayes object
varBayes = None
loaded_file_path = None

def load_varbayes_object(file_path):
  """Load VarBayes object from file"""
  global varBayes, loaded_file_path

  if not os.path.exists(file_path):
      raise FileNotFoundError(f"File not found: {file_path}")

  logger.info(f"Loading VarBayes object from: {file_path}")

  try:
      if file_path.endswith('.joblib'):
          varBayes = joblib.load(file_path)
      elif file_path.endswith(('.pickle', '.pkl')):
          with open(file_path, 'rb') as f:
              varBayes = pickle.load(f)
      else:
          raise ValueError("Unsupported file format. Use .joblib, .pickle, or .pkl")

      # Validate the object has the required methods
      if not hasattr(varBayes, 'check_cell') or not hasattr(varBayes, 'check_spot'):
          raise AttributeError("Object doesn't have required diagnostic methods (check_cell, check_spot)")

      loaded_file_path = file_path
      logger.info("✅ VarBayes object loaded successfully!")

  except Exception as e:
      logger.error(f"❌ Error loading VarBayes object: {e}")
      raise

@app.route('/check_cell', methods=['POST'])
def check_cell_endpoint():
  """Run check_cell diagnostic"""
  if varBayes is None:
      return jsonify({'success': False, 'error': 'No VarBayes object loaded'}), 500

  try:
      data = request.json
      if not data:
          return jsonify({'success': False, 'error': 'No JSON data provided'}), 400

      # Validate required parameters
      if 'cell_id' not in data or 'user_class' not in data:
          return jsonify({
              'success': False,
              'error': 'Missing required parameters: cell_id, user_class'
          }), 400

      cell_id = data['cell_id']
      user_class = data['user_class']
      top_n = data.get('top_n', 10)

      logger.info(f"Running check_cell for cell_id={cell_id}, user_class={user_class}")

      # Run diagnostic
      gene_data, fig = varBayes.check_cell(cell_id, user_class, top_n, show_plot=True)

      # Convert plot to base64
      buffer = io.BytesIO()
      fig.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
      buffer.seek(0)
      plot_data = base64.b64encode(buffer.read()).decode()
      plt.close(fig)

      return jsonify({
          'success': True,
          'plot_image': f"data:image/png;base64,{plot_data}",
          'gene_data': gene_data.fillna(None).to_dict('records'),
          'cell_id': cell_id,
          'user_class': user_class
      })

  except Exception as e:
      logger.error(f"Error in check_cell: {e}")
      return jsonify({'success': False, 'error': str(e)}), 500

@app.route('/check_spot', methods=['POST'])
def check_spot_endpoint():
    """Run check_spot diagnostic"""
    if varBayes is None:
        return jsonify({'success': False, 'error': 'No VarBayes object loaded'}), 500

    try:
        data = request.json
        logger.info(f"Received data: {data}")

        if not data:
            return jsonify({'success': False, 'error': 'No JSON data provided'}), 400

        if 'spot_id' not in data:
            return jsonify({
                'success': False,
                'error': 'Missing required parameter: spot_id'
            }), 400

        spot_id = data['spot_id']
        logger.info(f"Running check_spot for spot_id={spot_id}")

        # Run diagnostic
        logger.info("Calling varBayes.check_spot...")
        spot_data = varBayes.check_spot(spot_id)
        logger.info(f"Got spot_data type: {type(spot_data)}")
        logger.info(f"Spot_data shape: {spot_data.shape if hasattr(spot_data, 'shape') else 'no shape'}")

        # Handle NaN values by replacing them during dict conversion
        logger.info("Converting to dict and handling NaN values...")
        result_dict = spot_data.to_dict('records')

        # Replace NaN values with None in the resulting dictionaries
        import math
        for record in result_dict:
            for key, value in record.items():
                if isinstance(value, float) and math.isnan(value):
                    record[key] = None

        logger.info("Returning JSON...")
        return jsonify({
            'success': True,
            'spot_data': result_dict,
            'spot_id': spot_id
        })

    except Exception as e:
        logger.error(f"Error in check_spot: {e}")
        logger.error(f"Error type: {type(e)}")
        import traceback
        logger.error(f"Traceback: {traceback.format_exc()}")
        return jsonify({'success': False, 'error': str(e)}), 500


@app.route('/status', methods=['GET'])
def status():
  """Check server status"""
  return jsonify({
      'success': True,
      'varbayes_loaded': varBayes is not None,
      'loaded_file': loaded_file_path,
      'endpoints': ['check_cell', 'check_spot', 'status']
  })

def diagnostic_server(varbayes_path, host='0.0.0.0', port=5001, debug=False):
  """
  Start the pciSeq diagnostic server

  Args:
      varbayes_path (str): Path to VarBayes .joblib or .pickle file
      host (str): Host to bind to (default: '0.0.0.0')
      port (int): Port to bind to (default: 5001)
      debug (bool): Enable Flask debug mode (default: False)

  Example:
      import pciSeq
      pciSeq.diagnostic_server('/tmp/pciSeq/data/debug/pciSeq.joblib')
  """
  try:
      # Load the VarBayes object
      load_varbayes_object(varbayes_path)

      print(f"🚀 Starting pciSeq diagnostic server")
      print(f"📁 Loaded: {varbayes_path}")
      print(f"🌐 Server: http://{host}:{port}")
      print(f"📊 Available endpoints:")
      print(f"   POST http://{host}:{port}/check_cell")
      print(f"   POST http://{host}:{port}/check_spot")
      print(f"   GET  http://{host}:{port}/status")
      print(f"\n⏳ Server is running. Press Ctrl+C to stop.")

      # Start the Flask app
      app.run(host=host, port=port, debug=debug)

  except KeyboardInterrupt:
      print("\n👋 Server stopped by user")
  except Exception as e:
      logger.error(f"Failed to start server: {e}")
      raise


