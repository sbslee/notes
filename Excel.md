# Excel

## Frequently used commands <a name="Frequently-used-commands"></a>

To combine text from two cells:

If the cells `A1` and `B1` are `Hello` and `World`, respectively, then the command `=A1&" "&B1` will give `Hello World`.

To calculate current age from date of birth:

If the cell `A1` is a date in the format of MM/DD/YYYY, then the command `=INT((TODAY()-A1)/365)` will give the current age. For example, as of January 5, 2021, the command will output `33` for `10/9/1987`.
