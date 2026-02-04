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
            state = self._read_json(self.state_file)
            if "workflow_log" in state:
                del state["workflow_log"]
                self._write_json(self.state_file, state)
        else:
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
        try:
            with open(filepath, 'r') as f:
                return json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            return {}

    def _write_json(self, filepath, data):
        with open(filepath, 'w') as f:
            json.dump(data, f, indent=4)

    def log(self, message):
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        with open(self.log_file, 'a') as f:
            f.write(f"[{timestamp}] {message}\n")

    def update_state(self, category, data_dict):
        state = self._read_json(self.state_file)
        if category not in state:
            state[category] = {}
        state[category].update(data_dict)
        if "workflow_log" in state: del state["workflow_log"]
        self._write_json(self.state_file, state)
        self.log(f"Updated State: [{category}]")
        # print(f"📝 State updated: [{category}]")

    def save_result(self, metric, value, units=None):
        if metric == "simulation_status": return
        results = self._read_json(self.results_file)
        
        entry = value if isinstance(value, dict) else {"value": value}
        if units and not isinstance(value, dict): entry["units"] = units
            
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

    # --- ADDED THIS MISSING METHOD ---
    def get_result(self, metric):
        """Retrieves the value of a specific metric from results.json."""
        results = self._read_json(self.results_file)
        data = results.get(metric)
        
        # Handle format {"value": X, "units": Y}
        if isinstance(data, dict) and "value" in data:
            return data["value"]
        return data
    
    # Helper to access state properties safely
    @property
    def state(self):
        return self._read_json(self.state_file)
