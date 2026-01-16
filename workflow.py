import os
import json
from datetime import datetime

class WorkflowManager:
    """
    Manages state.json, results.json, and workflow.log.
    """
    def __init__(self, polymer_name):
        self.base_dir = os.path.abspath(polymer_name)
        os.makedirs(self.base_dir, exist_ok=True)
        
        self.state_file = os.path.join(self.base_dir, "state.json")
        self.results_file = os.path.join(self.base_dir, "results.json")
        self.log_file = os.path.join(self.base_dir, "workflow.log")
        
        # 1. Initialize or Clean State
        if os.path.exists(self.state_file):
            # If state exists, check for and REMOVE the old 'workflow_log' key
            # This fixes the issue where old logs persisted in the json.
            state = self._read_json(self.state_file)
            if "workflow_log" in state:
                del state["workflow_log"]
                self._write_json(self.state_file, state)
        else:
            # Create new clean state
            self._write_json(self.state_file, {
                "polymer_name": polymer_name,
                "created_at": datetime.now().isoformat(),
                "paths": {"base": self.base_dir}
            })
            
        # 2. Initialize Results
        if not os.path.exists(self.results_file):
            self._write_json(self.results_file, {})

        self.log("Workflow initialized.")

    def _read_json(self, filepath):
        with open(filepath, 'r') as f:
            return json.load(f)

    def _write_json(self, filepath, data):
        with open(filepath, 'w') as f:
            json.dump(data, f, indent=4)

    def log(self, message):
        """Appends a timestamped message to the external workflow.log file only."""
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        with open(self.log_file, 'a') as f:
            f.write(f"[{timestamp}] {message}\n")

    def update_state(self, category, data_dict):
        """Updates state.json without logging inside the json."""
        state = self._read_json(self.state_file)
        
        if category not in state:
            state[category] = {}
        
        state[category].update(data_dict)
        
        # Ensure we never accidentally write the log back
        if "workflow_log" in state:
            del state["workflow_log"]
            
        self._write_json(self.state_file, state)
        
        # Log to text file
        self.log(f"Updated State: [{category}]")
        print(f"📝 State updated: [{category}]")

    def save_result(self, metric, value, units=None):
        """
        Updates results.json. 
        Ensures 'degree_of_polymerization' is always the first entry.
        """
        if metric == "simulation_status":
            return

        results = self._read_json(self.results_file)
        
        # Create Entry
        if isinstance(value, dict):
             entry = value 
        else:
            entry = {"value": value}
            if units: entry["units"] = units
            
        results[metric] = entry
        
        # Re-order
        ordered_results = {}
        if "degree_of_polymerization" in results:
            ordered_results["degree_of_polymerization"] = results.pop("degree_of_polymerization")
        ordered_results.update(results)
        
        self._write_json(self.results_file, ordered_results)
        self.log(f"Result Saved: {metric}")
        print(f"📊 Result saved: {metric}")

    def get_path(self, key):
        state = self._read_json(self.state_file)
        return state.get("paths", {}).get(key)
    
    def add_path(self, key, path):
        self.update_state("paths", {key: os.path.abspath(path)})
