// ignore_for_file: deprecated_member_use, library_private_types_in_public_api

import 'package:flutter/material.dart';
import 'package:confetti/confetti.dart';
import 'package:google_fonts/google_fonts.dart';
import 'dart:convert';
import 'package:http/http.dart' as http;
import 'package:url_launcher/url_launcher.dart';
import 'package:flutter/services.dart';
import 'package:flutter_dotenv/flutter_dotenv.dart';

final String pythonApiUrl =
    dotenv.env['python_api_url'] ?? 'KEY_NOT_FOUND'; // Python API URL
final String geminiApiUrl = dotenv.env['geminiApiUrl'] ?? 'KEY_NOT_FOUND';
final String geminiApiKey =
    dotenv.env['GEMINI_api'] ?? 'KEY_NOT_FOUND'; // Replace with your API key
void main() async {
  await dotenv.load(fileName: '.env');
  runApp(const ChemMixApp());
}

class ChemMixApp extends StatelessWidget {
  const ChemMixApp({super.key});

  @override
  Widget build(BuildContext context) {
    return MaterialApp(
      title: 'ChemPredict Pro',
      debugShowCheckedModeBanner: false,
      theme: ThemeData(
        primaryColor: const Color(0xFF0A192F), // Navy blue
        colorScheme: ColorScheme.light(
          primary: const Color(0xFF0A192F),
          secondary: const Color(0xFF64FFDA), // Teal
          surface: const Color(0xFF172A45),
        ),
        scaffoldBackgroundColor: const Color(0xFF0A192F),
        textTheme: GoogleFonts.robotoTextTheme().copyWith(
          displayLarge: const TextStyle(
            color: Colors.white,
            fontWeight: FontWeight.bold,
          ),
          bodyLarge: const TextStyle(color: Colors.white70),
          bodyMedium: const TextStyle(color: Colors.white70),
        ),
        cardTheme: CardTheme(
          elevation: 8,
          shape: RoundedRectangleBorder(
            borderRadius: BorderRadius.circular(16),
          ),
          color: const Color(0xFF172A45),
        ),
        inputDecorationTheme: InputDecorationTheme(
          filled: true,
          fillColor: const Color(0xFF172A45).withOpacity(0.5),
          border: OutlineInputBorder(
            borderRadius: BorderRadius.circular(12),
            borderSide: BorderSide.none,
          ),
          contentPadding: const EdgeInsets.symmetric(
            horizontal: 20,
            vertical: 18,
          ),
        ),
      ),
      home: const HomePage(),
    );
  }
}

class HomePage extends StatefulWidget {
  const HomePage({super.key});

  @override
  _HomePageState createState() => _HomePageState();
}

class _HomePageState extends State<HomePage>
    with SingleTickerProviderStateMixin {
  late TabController _tabController;
  final TextEditingController smilesController = TextEditingController();
  final TextEditingController nameController = TextEditingController();
  bool isPredicting = false;
  bool showResult = false;
  late ConfettiController _confettiController;

  // Prediction results
  Map<String, dynamic> predictionResults = {};
  String geminiSummary = "";

  @override
  void initState() {
    super.initState();
    _tabController = TabController(length: 3, vsync: this);
    _confettiController =
        ConfettiController(duration: const Duration(seconds: 5));
  }

  @override
  void dispose() {
    _tabController.dispose();
    smilesController.dispose();
    nameController.dispose();
    _confettiController.dispose();
    super.dispose();
  }

  Future<Map<String, dynamic>> callPythonAPI(
      String smiles, String compoundName) async {
    try {
      Map<String, dynamic> data = {"smiles": smiles, "name": compoundName};
      final response = await http.post(
        Uri.parse(pythonApiUrl),
        headers: {"Content-Type": "application/json"},
        body: jsonEncode(data),
      );

      if (response.statusCode == 200) {
        return jsonDecode(response.body);
      } else {
        throw Exception(
            "Failed to fetch Python API results: ${response.statusCode}");
      }
    } catch (e) {
      throw Exception("Error calling Python API: $e");
    }
  }

  Future<String> callGeminiAPI(Map<String, dynamic> pythonResult) async {
    try {
      String prompt = buildGeminiPrompt(pythonResult);
      Map<String, dynamic> data = {
        "contents": [
          {
            "parts": [
              {"text": prompt}
            ]
          }
        ]
      };

      final response = await http.post(
        Uri.parse("$geminiApiUrl?key=$geminiApiKey"),
        headers: {"Content-Type": "application/json"},
        body: jsonEncode(data),
      );

      if (response.statusCode == 200) {
        Map<String, dynamic> result = jsonDecode(response.body);
        return result['candidates'][0]['content']['parts'][0]['text'];
      } else {
        throw Exception(
            "Failed to fetch Gemini API results: ${response.statusCode}");
      }
    } catch (e) {
      throw Exception("Error calling Gemini API: $e");
    }
  }

  String buildGeminiPrompt(Map<String, dynamic> pythonResult) {
    return """
    Analyze the following drug candidate:
    and give output in Four Lines without too technical ,all we need to know is if it is toxic are not
    Name: ${pythonResult['Name']}
    SMILES: ${pythonResult['SMILES']}
    Predicted Properties:
      - Molecular Weight: ${pythonResult['Predicted Properties']['MolWt']}
      - LogP: ${pythonResult['Predicted Properties']['MolLogP']}
      - TPSA: ${pythonResult['Predicted Properties']['TPSA']}
      - Hydrogen Donors: ${pythonResult['Predicted Properties']['HDonors']}
      - Hydrogen Acceptors: ${pythonResult['Predicted Properties']['HAcceptors']}
      - Ring Count: ${pythonResult['Predicted Properties']['RingCount']}
      - Rotatable Bonds: ${pythonResult['Predicted Properties']['Rotatable Bonds']}
      - Fraction CSP3: ${pythonResult['Predicted Properties']['FractionCSP3']}
    Toxicity Prediction: ${pythonResult['Toxicity Prediction']}
    Toxicity Confidence: ${pythonResult['Toxicity Confidence']}
    Provide a summary of the drug's safety and potential use , Please Dont add Further testing is needed.
    """;
  }

  Future<void> predict() async {
    if (smilesController.text.trim().isEmpty ||
        nameController.text.trim().isEmpty) {
      ScaffoldMessenger.of(context).showSnackBar(
        SnackBar(
          content: const Text(
              'Please fill in both SMILES notation and chemical name'),
          backgroundColor: Colors.red[800],
          behavior: SnackBarBehavior.floating,
          shape: RoundedRectangleBorder(
            borderRadius: BorderRadius.circular(12),
          ),
        ),
      );
      return;
    }

    setState(() {
      isPredicting = true;
      showResult = false;
    });

    try {
      // Step 1: Call Python API
      Map<String, dynamic> pythonResult = await callPythonAPI(
        smilesController.text,
        nameController.text,
      );
      predictionResults = pythonResult;

      // Step 2: Call Gemini API
      String geminiResult = await callGeminiAPI(pythonResult);
      geminiSummary = geminiResult;

      setState(() {
        showResult = true;
      });

      // Animate to results tab and play confetti
      _tabController.animateTo(1);
      _confettiController.play();
    } catch (e) {
      ScaffoldMessenger.of(context).showSnackBar(
        SnackBar(
          content: Text("Error: $e"),
          backgroundColor: Colors.red[800],
          behavior: SnackBarBehavior.floating,
          shape: RoundedRectangleBorder(
            borderRadius: BorderRadius.circular(12),
          ),
        ),
      );
    } finally {
      setState(() {
        isPredicting = false;
      });
    }
  }

  @override
  Widget build(BuildContext context) {
    return Stack(
      children: [
        Scaffold(
          backgroundColor: Colors.transparent,
          appBar: AppBar(
            title: Row(
              children: [
                const Icon(Icons.science, size: 30, color: Color(0xFF64FFDA)),
                const SizedBox(width: 12),
                Text(
                  'Drug Analyser',
                  style: GoogleFonts.orbitron(
                    fontSize: 24,
                    fontWeight: FontWeight.w700,
                    letterSpacing: 1.5,
                  ),
                ),
              ],
            ),
            centerTitle: false,
            elevation: 0,
            backgroundColor: Colors.transparent,
            bottom: TabBar(
              controller: _tabController,
              indicatorColor: const Color(0xFF64FFDA),
              indicatorWeight: 3,
              labelStyle: GoogleFonts.roboto(
                fontWeight: FontWeight.w500,
                fontSize: 16,
              ),
              unselectedLabelColor: Colors.white70,
              tabs: const [
                Tab(icon: Icon(Icons.science), text: 'Predict'),
                Tab(icon: Icon(Icons.analytics), text: 'Results'),
                Tab(icon: Icon(Icons.people), text: 'Team'),
              ],
            ),
          ),
          body: TabBarView(
            controller: _tabController,
            children: [
              buildPredictTab(),
              buildResultsTab(),
              buildTeamTab(),
            ],
          ),
        ),
        Align(
          alignment: Alignment.topCenter,
          child: ConfettiWidget(
            confettiController: _confettiController,
            blastDirectionality: BlastDirectionality.explosive,
            shouldLoop: false,
            colors: const [
              Color(0xFF64FFDA),
              Colors.white,
              Color(0xFF0A192F),
            ],
          ),
        ),
      ],
    );
  }

  Widget buildPredictTab() {
    return SingleChildScrollView(
      padding: const EdgeInsets.all(24),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          // Title
          Text(
            'Drug Analyser Engine',
            style: GoogleFonts.orbitron(
              fontSize: 28,
              fontWeight: FontWeight.bold,
              color: const Color(0xFF64FFDA),
            ),
          ),
          const SizedBox(height: 8),

          // Subtitle
          Text(
            'Enter chemical details to predict properties and interactions',
            style: GoogleFonts.roboto(
              fontSize: 16,
              color: Colors.white70,
            ),
          ),
          const SizedBox(height: 32),

          // Chemical Input Card
          Card(
            color: const Color(0xFF112240),
            elevation: 12,
            shape: RoundedRectangleBorder(
              borderRadius: BorderRadius.circular(16),
            ),
            child: Padding(
              padding: const EdgeInsets.all(24),
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  // SMILES Input
                  TextField(
                    controller: smilesController,
                    style: const TextStyle(color: Colors.white),
                    decoration: const InputDecoration(
                      labelText: 'SMILES Notation',
                      labelStyle: TextStyle(color: Colors.white70),
                      prefixIcon: Icon(Icons.code, color: Color(0xFF64FFDA)),
                      enabledBorder: UnderlineInputBorder(
                        borderSide: BorderSide(color: Color(0xFF64FFDA)),
                      ),
                    ),
                  ),
                  const SizedBox(height: 24),

                  // Chemical Name Input
                  TextField(
                    controller: nameController,
                    style: const TextStyle(color: Colors.white),
                    decoration: const InputDecoration(
                      labelText: 'Chemical Name',
                      labelStyle: TextStyle(color: Colors.white70),
                      prefixIcon:
                          Icon(Icons.emoji_objects, color: Color(0xFF64FFDA)),
                      enabledBorder: UnderlineInputBorder(
                        borderSide: BorderSide(color: Color(0xFF64FFDA)),
                      ),
                    ),
                  ),
                  const SizedBox(height: 32),

                  // Predict Button
                  SizedBox(
                    width: double.infinity,
                    child: ElevatedButton.icon(
                      onPressed: isPredicting ? null : predict,
                      icon: const Icon(Icons.science),
                      label: isPredicting
                          ? const CircularProgressIndicator(
                              color: Colors.white,
                            )
                          : const Text(
                              "Predict Properties",
                              style: TextStyle(fontSize: 16),
                            ),
                      style: ElevatedButton.styleFrom(
                        backgroundColor: const Color(0xFF64FFDA),
                        foregroundColor: Colors.black,
                        padding: const EdgeInsets.symmetric(vertical: 16),
                        shape: RoundedRectangleBorder(
                          borderRadius: BorderRadius.circular(12),
                        ),
                      ),
                    ),
                  ),
                ],
              ),
            ),
          ),

          const SizedBox(height: 32),

// 💡Example card
          Card(
            color: const Color(0xFF0F172A),
            shape: RoundedRectangleBorder(
              borderRadius: BorderRadius.circular(16),
            ),
            elevation: 6,
            child: Padding(
              padding: const EdgeInsets.all(20),
              child: Row(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  const Icon(Icons.lightbulb,
                      color: Color(0xFF64FFDA), size: 30),
                  const SizedBox(width: 16),
                  Expanded(
                    child: Text(
                      "EXAMPLE INPUTS: SMILES Notation : C1=CC=C2C(=C1)C(=O)N(C2=O)C3=CC=CC=C3 , CHEMICAL NAME :CHEMBL1234567",
                      style: GoogleFonts.robotoMono(
                          fontSize: 10, color: Colors.white70),
                    ),
                  ),
                ],
              ),
            ),
          ),
          const SizedBox(height: 32),

// 💡 Fun Chemistry Fact 1
          Card(
            color: const Color(0xFF0F172A),
            shape: RoundedRectangleBorder(
              borderRadius: BorderRadius.circular(16),
            ),
            elevation: 6,
            child: Padding(
              padding: const EdgeInsets.all(20),
              child: Row(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  const Icon(Icons.lightbulb,
                      color: Color(0xFF64FFDA), size: 30),
                  const SizedBox(width: 16),
                  Expanded(
                    child: Text(
                      "💡 Did you know? Mixing hydrogen peroxide with potassium iodide triggers a spectacular reaction known as the 'Elephant's Toothpaste' – releasing rapid foamy oxygen!",
                      style: GoogleFonts.robotoMono(
                          fontSize: 14, color: Colors.white70),
                    ),
                  ),
                ],
              ),
            ),
          ),

          const SizedBox(height: 16),

// 💡 Fact 3: Lipinski's Rule
          Card(
            color: const Color(0xFF0F172A),
            shape: RoundedRectangleBorder(
              borderRadius: BorderRadius.circular(16),
            ),
            elevation: 6,
            child: Padding(
              padding: const EdgeInsets.all(20),
              child: Row(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  const Icon(Icons.rule, color: Color(0xFF64FFDA), size: 30),
                  const SizedBox(width: 16),
                  Expanded(
                    child: Text(
                      "📏 Did You know ? Lipinski’s Rule of Five is a rule of thumb to evaluate drug-likeness. It predicts good oral bioavailability if the compound has ≤5 H-bond donors, ≤10 H-bond acceptors, MW < 500, and logP < 5.",
                      style: GoogleFonts.robotoMono(
                          fontSize: 14, color: Colors.white70),
                    ),
                  ),
                ],
              ),
            ),
          ),

          const SizedBox(height: 16),

// 💡 Fact 4: Molecular Docking
          Card(
            color: const Color(0xFF0F172A),
            shape: RoundedRectangleBorder(
              borderRadius: BorderRadius.circular(16),
            ),
            elevation: 6,
            child: Padding(
              padding: const EdgeInsets.all(20),
              child: Row(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  const Icon(Icons.biotech, color: Color(0xFF64FFDA), size: 30),
                  const SizedBox(width: 16),
                  Expanded(
                    child: Text(
                      "⚙️Did You know? Molecular docking predicts how two molecules, such as a drug and protein, might bind together. It’s a vital tool in computational drug design for estimating binding affinity and stability.",
                      style: GoogleFonts.robotoMono(
                          fontSize: 14, color: Colors.white70),
                    ),
                  ),
                ],
              ),
            ),
          ),
        ],
      ),
    );
  }

  String formatPredictionResults(Map<String, dynamic> results) {
    final buffer = StringBuffer();

    buffer.writeln('Name: ${results['Name']}');

    final props = results['Predicted Properties'] as Map<String, dynamic>;
    buffer.writeln('--- Predicted Properties ---');
    props.forEach((key, value) {
      buffer.writeln('$key: $value');
    });

    buffer.writeln('SMILES: ${results['SMILES']}');
    buffer.writeln('Toxicity Confidence: ${results['Toxicity Confidence']}');
    buffer.writeln('Toxicity Prediction: ${results['Toxicity Prediction']}');

    return buffer.toString();
  }

  Widget buildResultsTab() {
    if (!showResult) {
      return Center(
        child: Column(
          mainAxisAlignment: MainAxisAlignment.center,
          children: [
            Icon(
              Icons.science,
              size: 80,
              color: Colors.white.withOpacity(0.3),
            ),
            const SizedBox(height: 24),
            Text(
              'Enter chemical details to see predictions',
              style: GoogleFonts.roboto(
                fontSize: 18,
                color: Colors.white70,
              ),
            ),
          ],
        ),
      );
    }

    return SingleChildScrollView(
      padding: const EdgeInsets.all(24),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          Text(
            'Prediction Results',
            style: GoogleFonts.orbitron(
              fontSize: 28,
              fontWeight: FontWeight.w700,
              color: const Color(0xFF64FFDA),
            ),
          ),
          const SizedBox(height: 16),
          Text(
            formatPredictionResults(predictionResults),
            style: GoogleFonts.robotoMono(fontSize: 20, color: Colors.white70),
          ),
          const SizedBox(height: 32),
          Text(
            'Gemini Summary',
            style: GoogleFonts.orbitron(
              fontSize: 30,
              fontWeight: FontWeight.w700,
              color: const Color(0xFF64FFDA),
            ),
          ),
          const SizedBox(height: 16),
          Text(
            geminiSummary,
            style: GoogleFonts.robotoMono(fontSize: 20, color: Colors.white70),
          ),
        ],
      ),
    );
  }

  Widget buildTeamTab() {
    return SingleChildScrollView(
      padding: const EdgeInsets.all(24),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          Text(
            'Our Details',
            style: GoogleFonts.orbitron(
              fontSize: 32,
              fontWeight: FontWeight.w700,
              color: const Color(0xFF64FFDA),
            ),
          ),
          const SizedBox(height: 16),
          const Divider(color: Color(0xFF64FFDA), thickness: 1),
          const SizedBox(height: 24),
          _buildProjectLinks(),
          const SizedBox(height: 40),
          Text(
            'Team Members',
            style: GoogleFonts.orbitron(
              fontSize: 24,
              fontWeight: FontWeight.w600,
              color: Colors.white,
            ),
          ),
          const SizedBox(height: 24),
          GridView.count(
            shrinkWrap: true,
            physics: const NeverScrollableScrollPhysics(),
            crossAxisCount: 2,
            crossAxisSpacing: 24,
            mainAxisSpacing: 24,
            childAspectRatio: 1.0,
            children: [
              _buildTeamMember(
                  name: 'Ajay Krishna',
                  role: 'ML Engineer',
                  github: 'https://github.com/ajaykrishna00-7',
                  linkedin:
                      'https://www.linkedin.com/in/ajay-krishna-588a1228b/',
                  imagePath: 'assets/Profile/ajay.jpg'),
              _buildTeamMember(
                  name: 'Dashrad Raghav',
                  role: 'ML Engineer',
                  github: 'https://github.com/cheetahraghav',
                  linkedin: 'https://www.linkedin.com/in/dhasradraghavab/',
                  imagePath: 'assets/Profile/raghav.jpg'),
              _buildTeamMember(
                  name: 'Bargavan',
                  role: 'Frontend/Backend Engineer',
                  github: 'https://github.com/BargavanR',
                  linkedin: 'https://www.linkedin.com/in/bargavan/',
                  imagePath: 'assets/Profile/bargavan.jpg'),
              _buildTeamMember(
                  name: 'Arun',
                  role: 'Frontend Engineer',
                  github: '',
                  linkedin: 'https://www.linkedin.com/in/arun-r-012b39277/',
                  imagePath: 'assets/Profile/arun.jpeg'),
            ],
          ),
        ],
      ),
    );
  }

  Widget _buildProjectLinks() {
    return Container(
      padding: const EdgeInsets.all(20),
      decoration: BoxDecoration(
        color: const Color(0xFF172A45),
        borderRadius: BorderRadius.circular(16),
        boxShadow: [
          BoxShadow(
            color: Colors.black26,
            blurRadius: 10,
            offset: const Offset(0, 4),
          ),
        ],
      ),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          Text(
            'Project Resources',
            style: GoogleFonts.roboto(
              fontSize: 20,
              fontWeight: FontWeight.bold,
              color: Colors.white,
            ),
          ),
          const SizedBox(height: 16),
          _buildResourceLink(
            icon: Icons.description,
            title: 'Project Presentation',
            subtitle: 'DrugAnalyser Pro concept and implementation details',
            url:
                'https://annauniv0-my.sharepoint.com/:f:/g/personal/2023511001_student_annauniv_edu/Eq0ay1FuL0tHnzhhBv8AFn4Bfy8YuhTn4Qz00QS8kEKFyA?e=K4LHeO',
          ),
          const Divider(color: Colors.white24, height: 24),
          _buildResourceLink(
            icon: Icons.code,
            title: 'GitHub Repository',
            subtitle: 'Source code and documentation',
            url: 'https://github.com/BargavanR/Drug_Analyser-GDG-2025',
          ),
        ],
      ),
    );
  }

  void openUrl(String url) async {
    final Uri uri = Uri.parse(url);
    if (await canLaunchUrl(uri)) {
      await launchUrl(uri, mode: LaunchMode.externalApplication);
    } else {
      throw 'Could not launch $url';
    }
  }

  Widget _buildResourceLink({
    required IconData icon,
    required String title,
    required String subtitle,
    required String url,
  }) {
    return InkWell(
      onTap: () => launchUrl(Uri.parse(url)),
      child: Padding(
        padding: const EdgeInsets.symmetric(vertical: 8),
        child: Row(
          children: [
            Container(
              width: 48,
              height: 48,
              decoration: BoxDecoration(
                color: const Color(0xFF0A192F),
                borderRadius: BorderRadius.circular(12),
              ),
              child: Icon(icon, color: const Color(0xFF64FFDA), size: 24),
            ),
            const SizedBox(width: 16),
            Expanded(
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  Text(
                    title,
                    style: GoogleFonts.roboto(
                      fontSize: 20,
                      fontWeight: FontWeight.bold,
                      color: Colors.white,
                    ),
                  ),
                  Text(
                    subtitle,
                    style: GoogleFonts.roboto(
                      fontSize: 15,
                      color: Colors.white70,
                    ),
                  ),
                ],
              ),
            ),
            Icon(Icons.open_in_new, color: const Color(0xFF64FFDA), size: 20),
          ],
        ),
      ),
    );
  }

  Widget _buildTeamMember({
    required String name,
    required String role,
    required String github,
    required String linkedin,
    required String imagePath, // <-- New parameter
  }) {
    return Card(
      elevation: 10,
      shape: RoundedRectangleBorder(
        borderRadius: BorderRadius.circular(16),
      ),
      color: const Color(0xFF0A192F),
      child: Padding(
        padding: const EdgeInsets.all(20),
        child: Column(
          mainAxisAlignment: MainAxisAlignment.center,
          children: [
            Container(
              width: 400,
              height: 400,
              decoration: BoxDecoration(
                shape: BoxShape.circle,
                border: Border.all(
                  color: const Color(0xFF64FFDA),
                  width: 2,
                ),
                color: const Color(0xFF172A45),
              ),
              child: ClipOval(
                child: Image.asset(
                  imagePath,
                  fit: BoxFit.cover,
                ),
              ),
            ),
            const SizedBox(height: 16),
            Text(
              name,
              style: GoogleFonts.roboto(
                fontSize: 18,
                fontWeight: FontWeight.bold,
                color: Colors.white,
              ),
            ),
            const SizedBox(height: 4),
            Text(
              role,
              style: GoogleFonts.roboto(
                fontSize: 14,
                color: const Color(0xFF64FFDA),
              ),
            ),
            const SizedBox(height: 16),
            Row(
              mainAxisAlignment: MainAxisAlignment.center,
              children: [
                IconButton(
                  icon: const Icon(Icons.link),
                  color: Colors.white,
                  onPressed: () => openUrl(linkedin),
                  tooltip: 'LinkedIn Profile',
                ),
                IconButton(
                  icon: const Icon(Icons.code),
                  color: Colors.white,
                  onPressed: () => openUrl(github),
                  tooltip: 'GitHub Profile',
                ),
              ],
            ),
          ],
        ),
      ),
    );
  }
}
