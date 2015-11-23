/* Copyright (C) 2015 Erik Källman */
/* This file is part of std_char_ops. */

/* std_char_ops is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU Lesser General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* std_char_ops is distributed in the hope that it will be useful, */
/* but without any warranty; without even the implied warranty of */
/* merchantability or fitness for a particular purpose. See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with std_char_ops, found in the "license" subdirectory of the root */
/* directory of any program using the std_char_ops library.*/
/*   If not, see <http://www.gnu.org/licenses/>. */
/**
   * @file std_char_ops.h
   * @author Erik Källman
   * @date November 2015
   * @brief This file contains public interface to the functions in the
   * std_char_opt library.
   */
#ifndef CHAR_OPS_H
#define CHAR_OPS_H

/**
   * @brief The concs() function concatenates @p n_args number of strings into
   * one string. it is convenient for creating string for directory paths that
   * are made up of different numbers of substrings.
   *
   * @param n_args The number of strings to be concatenated.
   * @param ... The strings to be concatenated.
   * @returns path The string containing the concatenated strings.
   * @note Uses variadic arguments, example call:
   * char *new_string = concs(3,string1, string2, string3);
   * @note Side effects: allocates memory for the string @p path to be returned.
   * by the function.
   */
char *
concs (int n_args, ...);

/**
   * @brief Use the standard UNIX strstr function to check if a string @p str1
   * is present in any of the strings contained in @p str2.
   *
   * @param str1 The string to look for.
   * @param str2 The two-dimensional array of string in which to look for
   * @p str1.
   * @returns The index @p j of the element in @p str2 that contains the
   * substring str1. If j is not found, returns @p -1.
   */
int
strscmp (const char *str1, const char **str2, int n_str);

/**
   * @brief This function uses a reference key @p num_key to determine whether
   * or not any of the chars in the string @p fall into a certain category of
   * alphabetical and symbolic characters.
   *
   * @param s The string to check for the given character type.
   * @param len The length of the string @p s in characters.
   * @returns 1 if the character type is found, otherwise 0.
   */
int
isanyalpha (char *s, int len);

/**
   * @brief Checks if a certain character @p is contained in a string @p str.
   *
   * @param str The string in which to search for @p c.
   * @param c The character to search for.
   */
int
charinstr (char *str, char c);

/**
   * @brief Extract the power of a number contained in the string @p s.
   *
   * @param s The string containing the number.
   * @param len The length of the string @p s.
   * @returns The power of the number contained in s, as an integer value.
   * @note Exits with EXIT_FAILURE if the power of the number is larger than
   * @p MAX_POWERL
   */
int
satopow(char *s, int len);

/**
   * @brief Extract a number in scientific notation, contained in a string @p s,
   * to a double.
   *
   * This function is useful for cases when numbers in scientific notation are
   * read from file as strings, and need be convert for the program to use
   * the numbers in subsequent calculations.
   *
   * @param s The string containing the number to be extracted and converted.
   * @param len Length of the string @p s in characters.
   * @returns The value of the number inside @p s, converted to @p double type.
   */
double
satof(char *s, int len);

/**
   * @brief This function converts an array of characters to @p double type
   * number even if its in scientific notation.
   *
   * @param s The string containing the number to be extracted and converted.
   * @param len The length of the string @p s, in characters.
   * @returns The value of the scientifically notated number, extracted from @å s.
   */
double
sci_atof(char *s);

/**
   * @brief Check if the string @p s contains only dashes.
   *
   * @param s The string to check for dashes.
   * @returns 1 if at least one double dash is found, otherwise 0, if none
   are found.
   */
int
isdashes (char *s, int len);

/**
   * @brief Check if tghe string @p s contains nothing but "empty space"
   characters (tabs, newlines, blankspace, etc.).
   *
   * @param s The string to be checked if "empty".
   * @param len The length of the string @p s, in characters.
   * @returns 0 if a character not defined in @p key is found, otherwise
   * returns 1.
   */
int
isempty (char *s, int len);

/**
   * @brief Checks if the string @p s contains any digit characters.
   *
   * @param s The string to check for digit characters.
   * @param len The length of the string @p s, in characters.
   * @returns 1 if digit found, otherwise 0.
   */
int
isdigitin (char *s, int len);

#endif /* CHAR_OPS_H */
