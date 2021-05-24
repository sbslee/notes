Excel
*****

Frequently used commands for Excel
==================================

* To combine text from two cells:

    If the cells ``A1`` and ``B1`` are ``Hello`` and ``World``, respectively, then the command ``=A1&" "&B1`` will give ``Hello World``.

* To calculate current age from date of birth:

    If the cell ``A1`` is a date in the format of MM/DD/YYYY, then the command ``=INT((TODAY()-A1)/365)`` will give the current age. For example, as of January 5, 2021, the command will output ``33`` for ``10/9/1987``.

* To change the date format:

    If the cell ``A1`` is a date in the format of MM/DD/YYYY, then the command ``=TEXT(A1,"yyyy-mm")`` will give YYYY-DD.

* To count asterisks in a cell:

    If the cell ``A1`` is ``Hello**World***``, then the command ``=SUM(LEN(A1)-LEN(SUBSTITUTE(A1,"*","")))`` will give ``5``.

* To remove text before or after delimiter:

    If the cell ``A1`` is ``Hello World``, then the command ``=LEFT(A1,FIND(" ",A1)-1)`` will give ``Hello``. If you want ``World``, then use the command ``=RIGHT(A1,LEN(A1)-FIND(" ",A1))``.

* To highlight duplicates in Google Sheets:

    ``countif(A:A,A1)>1``

* To count cells with fill color in Excel (Developer > Macros):

    .. code-block:: console

				Function GetColor(r As Range) As Integer
					GetColor = r.Interior.ColorIndex
				End Function

				Function CountColor(MyRange As Range, ColorIndex As Integer)
					Dim iCount As Integer
					Application.Volatile
					iCount = 0
					For Each cell In MyRange
						If cell.Interior.ColorIndex = ColorIndex Then
							iCount = iCount + 1
						End If
					Next cell
					CountColor = iCount
				End Function

References:

  - `How to calculate age in Excel: from date of birth, between two dates, and more <https://www.ablebits.com/office-addins-blog/2016/10/19/calculate-age-excel/#:~:text=Simply%20by%20subtracting%20the%20birth,also%20be%20used%20in%20Excel.&text=The%20first%20part%20of%20the,get%20the%20numbers%20of%20years.>`__
  - `How To Highlight/Find And Remove Duplicates In Google Sheets <https://www.alphr.com/highlight-duplicates-google-sheets/>`__
