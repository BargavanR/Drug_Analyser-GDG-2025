// This is a basic Flutter widget test.
//
// To perform an interaction with a widget in your test, use the WidgetTester
// utility in the flutter_test package. For example, you can send tap and scroll
// gestures. You can also use WidgetTester to find child widgets in the widget
// tree, read text, and verify that the values of widget properties are correct.

// ignore_for_file: unused_import

import 'package:flutter/material.dart';
import 'package:flutter_test/flutter_test.dart';

import 'package:chemix88/main.dart';

void main() {
  testWidgets('Counter increments smoke test', (WidgetTester tester) async {
    // Build our app and trigger a frame.
    await tester.pumpWidget(const ChemMixApp());

    // Verify that the app title is present
    expect(find.text('ChemPredict Pro'), findsOneWidget);

    // Verify that the main tabs are present
    expect(find.text('Predict'), findsOneWidget);
    expect(find.text('Results'), findsOneWidget);
    expect(find.text('Analysis'), findsOneWidget);
    expect(find.text('Team'), findsOneWidget);
  });
}
