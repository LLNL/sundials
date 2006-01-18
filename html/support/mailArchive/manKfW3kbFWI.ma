(************** Content-type: application/mathematica **************

                    Mathematica-Compatible Notebook

This notebook can be used with any Mathematica-compatible
application, such as Mathematica, MathReader or Publicon. The data
for the notebook starts with the line containing stars above.

To get the notebook into a Mathematica-compatible application, do
one of the following:

* Save the data starting with the line of stars above into a file
  with a name ending in .nb, then open the file inside the
  application;

* Copy the data starting with the line of stars above to the
  clipboard, then use the Paste menu command inside the application.

Data for notebooks contains only printable 7-bit ASCII and can be
sent directly in email or through ftp in text mode.  Newlines can be
CR, LF or CRLF (Unix, Macintosh or MS-DOS style).

NOTE: If you modify the data for this notebook not in a Mathematica-
compatible application, you must delete the line below containing
the word CacheID, otherwise Mathematica-compatible applications may
try to use invalid cache data.

For more information on notebooks and Mathematica-compatible 
applications, contact Wolfram Research:
  web: http://www.wolfram.com
  email: info@wolfram.com
  phone: +1-217-398-0700 (U.S.)

Notebook reader applications are available free of charge from 
Wolfram Research.
*******************************************************************)

(*CacheID: 232*)


(*NotebookFileLineBreakTest
NotebookFileLineBreakTest*)
(*NotebookOptionsPosition[      3094,         83]*)
(*NotebookOutlinePosition[      3737,        105]*)
(*  CellTagsIndexPosition[      3693,        101]*)
(*WindowFrame->Normal*)



Notebook[{
Cell[BoxData[
    \(NDSolve[{\((x^2/4 + 2 \((1 - y[x]\ y[x])\))\) \(y''\)[x] + 
            x^2  y[x]\ \(y'\)[x]^2/\((4 \((1 - y[x]^2)\))\) + 
            x\ \(y'\)[x]/2 + 
            2  y[x] \((1 - 
                  y[x]\ y[x])\) \((1/4 + \((1 - y[x]^2)\)/\((x^2)\))\) == 0, 
        y[0] \[Equal] 0, y[\[Infinity]] \[Equal] 0}, 
      y[x], {x, 0, \[Infinity]}]\)], "Input"],

Cell[BoxData[{
    RowBox[{"\[IndentingNewLine]", 
      RowBox[{
        RowBox[{\(2\ y[
              x]\ \((1 - y[x]\^2)\)\ \((1\/4 + \(1 - y[x]\^2\)\/x\^2)\)\), 
          "+", 
          RowBox[{\(1\/2\), " ", "x", " ", 
            RowBox[{
              SuperscriptBox["y", "\[Prime]",
                MultilineFunction->None], "[", "x", "]"}]}], "+", 
          FractionBox[
            RowBox[{\(x\^2\), " ", \(y[x]\), " ", 
              SuperscriptBox[
                RowBox[{
                  SuperscriptBox["y", "\[Prime]",
                    MultilineFunction->None], "[", "x", "]"}], 
                "2"]}], \(4\ \((1 - y[x]\^2)\)\)], "+", 
          RowBox[{\((x\^2\/4 + 2\ \((1 - y[x]\^2)\))\), " ", 
            RowBox[{
              SuperscriptBox["y", "\[Prime]\[Prime]",
                MultilineFunction->None], "[", "x", "]"}]}]}], "\[Equal]", 
        "0"}]}], "\[IndentingNewLine]", \(y[0] \[Equal] 0, 
    y[\[Infinity]] \[Equal] 0\), "\[IndentingNewLine]", " "}], "Input"]
},
FrontEndVersion->"4.1 for Microsoft Windows",
ScreenRectangle->{{0, 1024}, {0, 681}},
WindowSize->{913, 523},
WindowMargins->{{0, Automatic}, {Automatic, 0}}
]

(*******************************************************************
Cached data follows.  If you edit this Notebook file directly, not
using Mathematica, you must remove the line containing CacheID at
the top of  the file.  The cache data will then be recreated when
you save this file from within Mathematica.
*******************************************************************)

(*CellTagsOutline
CellTagsIndex->{}
*)

(*CellTagsIndex
CellTagsIndex->{}
*)

(*NotebookFileOutline
Notebook[{
Cell[1705, 50, 378, 7, 70, "Input"],
Cell[2086, 59, 1004, 22, 110, "Input"]
}
]
*)



(*******************************************************************
End of Mathematica Notebook file.
*******************************************************************)

