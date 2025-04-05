from flask import Flask, request, jsonify
from rdkit import Chem
from rdkit.Chem import Descriptors
import joblib
import json

# Load pre-trained model and scaler
final_model = joblib.load(r"E:\Bargavan\Projects\GDG\chemix88_v1\chemix88\model_files\drug_toxicity_model.pkl")
scaler = joblib.load(r"E:\Bargavan\Projects\GDG\chemix88_v1\chemix88\model_files\scaler.pkl")

app = Flask(__name__)

def extract_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return [0] * 8
    return [
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        Descriptors.TPSA(mol),
        Descriptors.NumHDonors(mol),
        Descriptors.NumHAcceptors(mol),
        Descriptors.RingCount(mol),
        Descriptors.NumRotatableBonds(mol),
        Descriptors.FractionCSP3(mol),
    ]

@app.route('/predict', methods=['POST'])
def predict():
    data = request.json
    smiles = data.get("smiles")
    compound_name = data.get("name", "")

    desc = extract_descriptors(smiles)
    scaled = scaler.transform([desc])
    prob = final_model.predict_proba(scaled)[0]
    pred = final_model.predict(scaled)[0]

    result = {
        "Name": compound_name,
        "SMILES": smiles,
        "Predicted Properties": {
            "MolWt": desc[0],
            "MolLogP": desc[1],
            "TPSA": desc[2],
            "HDonors": desc[3],
            "HAcceptors": desc[4],
            "RingCount": desc[5],
            "Rotatable Bonds": desc[6],
            "FractionCSP3": desc[7],
        },
        "Toxicity Prediction": int(pred),
        "Toxicity Confidence": float(prob[1]),
    }
    return jsonify(result)

if __name__ == "__main__":
    app.run(debug=True)