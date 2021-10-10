Miscellaneous
*************

NET::ERR_CERT_INVALID
=====================

Google Chrome will throw this error if the website uses https which you haven’t added to your approved certs. This annoying error can be overcome if you type “badidea” or “thisisunsafe” directly in Chrome on the same page.

References:

  - `Chrome: Bypass NET::ERR_CERT_INVALID for development <https://medium.com/@dblazeski/chrome-bypass-net-err-cert-invalid-for-development-daefae43eb12>`__

Terminal
========

Moving the cursor
-----------------

+------------------+-------------------------------------------+
| Keys             | Cursor behavior                           |
+==================+===========================================+
| ``Esc`` + ``F``  | Move the cursor forward by one word       |
+------------------+-------------------------------------------+
| ``Esc`` + ``B``  | Move the cursor backward by one word      |
+------------------+-------------------------------------------+
| ``Ctrl`` + ``A`` | Jump to the beginning of the line         |
+------------------+-------------------------------------------+
| ``Ctrl`` + ``E`` | Jump to the end of the line               |
+------------------+-------------------------------------------+
| ``Ctrl`` + ``U`` | Clear the entire part left of the cursor  |
+------------------+-------------------------------------------+
| ``Ctrl`` + ``K`` | Clear the entire part right of the cursor |
+------------------+-------------------------------------------+
| ``Ctrl`` + ``C`` | Abort and go to the new line              |
+------------------+-------------------------------------------+


Google Foobar
=============

Lovely Lucky LAMBs
------------------

`Reference <https://datanonymous.wordpress.com/foobar-level-2-lovely-lucky-lambs/>`__

.. code-block:: note

    def solution(total_lambs):

        generous_list = []
        n = 0
        counter = 0

        while counter < 1E9:
            next = 2 ** n
            if counter + next > total_lambs:
                break
            counter += next
            generous_list.append(next)
            n += 1

        if total_lambs == 1:
            stingy_list = [1]
            n = 1
            counter = 1
        elif total_lambs == 2:
            stingy_list = [1, 1]
            n = 2
            counter = 2
        else:
            stingy_list = [1, 1]
            n = 2
            counter = 2

            while counter <= 1E9:
                next = stingy_list[n-1] + stingy_list[n-2]
                if counter + next > total_lambs:
                    break
                stingy_list.append(next)
                counter += next
                n += 1

        return len(stingy_list) - len(generous_list)

    print solution(10)  # Gives 1 - [1, 1, 2, 3] vs. [1, 2, 4]
    print solution(143) # Gives 3 - [1, 1, 2, 3, 5, 8, 13, 21, 34, 55] vs. [1, 2, 4, 8, 16, 32, 64]

Find the Access Codes
---------------------

`Reference <https://stackoverflow.com/questions/39846735/google-foobar-challenge-3-find-the-access-codes>`__

.. code-block:: note

    def solution(l):
        c = [0] * len(l)
        count = 0
        for i in range(len(l)):
            for j in range(i):
                if l[i] % l[j] == 0:
                    c[i] += 1
                    count += c[j]
        return count

    print solution([1, 1, 1])              # Gives 1 - [1, 1, 1]
    print solution([1, 2, 3, 4, 5, 6])     # Gives 3 - [1, 2, 4], [1, 2, 6], [1, 3, 6]
    print solution([1, 2, 3, 4, 5, 6, 10]) # Gives 5 - [1, 2, 4], [1, 2, 6], [1, 2, 10], [1, 3, 6], [1, 5, 10]
