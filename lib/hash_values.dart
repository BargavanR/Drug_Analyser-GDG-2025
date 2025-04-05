import 'package:flutter/foundation.dart';

// Provide a polyfill for the deprecated hashValues method
int hashValues(Object? a, [Object? b, Object? c, Object? d, Object? e]) {
  return Object.hash(a, b, c, d, e);
}
