let model;
let scaler;

async function loadModel() {
    model = await tf.loadLayersModel('assets/model.json');
    const response = await fetch('assets/scaler.json');
    scaler = await response.json();
}

async function predict(inputData) {
    // Scale the input data using the scaler
    const scaledData = inputData.map((value, index) => value / scaler[index]);
    const inputTensor = tf.tensor2d([scaledData]);
    const prediction = model.predict(inputTensor);
    return prediction.dataSync(); // Return prediction as an array
}