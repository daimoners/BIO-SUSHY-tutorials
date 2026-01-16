import os
import json
from datetime import datetime

class WorkflowManager:
    """
    Manages the state.json and results.json files for the experiment.
    """
    def __init__(self, polymer_name):
        self.base_dir = os.path.abspath(polymer_name)
        os.makedirs(self.base_dir, exist_ok=True)
        
        self.state_file = os.path.join(self.base_dir, "state.json")
        self.results_file = os.path.join(self.base_dir, "results.json")
        
        # Initialize if not exists
        if not os.path.exists(self.state_file):
            self._write_json(self.state_file, {
                "polymer_name": polymer_name,
                "created_at": datetime.now().isoformat(),
                "status": "initialized",
                "paths": {"base": self.base_dir},
                "workflow_log": []
            })
            
        if not os.path.exists(self.results_file):
            self._write_json(self.results_file, {})

    def _read_json(self, filepath):
        with open(filepath, 'r') as f:
            return json.load(f)

    def _write_json(self, filepath, data):
        with open(filepath, 'w') as f:
            json.dump(data, f, indent=4)

    def update_state(self, category, data_dict):
        """
        Updates state.json with metadata (paths, settings, status).
        Args:
            category (str): The section name (e.g., 'simulation_settings', 'build_step').
            data_dict (dict): The data to store.
        """
        state = self._read_json(self.state_file)
        
        # Ensure category exists or update it
        if category not in state:
            state[category] = {}
        
        # Merge dictionaries (shallow merge)
        state[category].update(data_dict)
        
        # Log update
        state["workflow_log"].append(f"{datetime.now().strftime('%H:%M:%S')} - Updated {category}")
        
        self._write_json(self.state_file, state)
        print(f"📝 State updated: [{category}]")

    def save_result(self, metric, value, units=None):
        """
        Updates results.json with scientific data (atom counts, energies, etc).
        """
        results = self._read_json(self.results_file)
        
        entry = {"value": value}
        if units:
            entry["units"] = units
            
        results[metric] = entry
        self._write_json(self.results_file, results)
        print(f"📊 Result saved: {metric} = {value} {units if units else ''}")

    def get_path(self, key):
        """Retrives a path stored in the state under 'paths'"""
        state = self._read_json(self.state_file)
        return state.get("paths", {}).get(key)
    
    def add_path(self, key, path):
        """Adds a file path to the 'paths' section of the state"""
        self.update_state("paths", {key: os.path.abspath(path)})
