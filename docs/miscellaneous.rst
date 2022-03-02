Miscellaneous
*************

scRNAseq analysis
=================

- `[Biostars] Single cell RNA-seq: how many PCs to use for t-SNE/UMAP? <https://www.biostars.org/p/9466164/>`__

Stranded vs Non-Stranded RNAseq
===============================

This `article <https://blog.genewiz.com/stranded-vs-non-stranded-rna-seq?utm_term=&utm_campaign=*APAC*+NGS&utm_source=adwords&utm_medium=ppc&hsa_tgt=dsa-460355903923&hsa_grp=93772107045&hsa_src=g&hsa_net=adwords&hsa_mt=&hsa_ver=3&hsa_ad=419449356856&hsa_acc=8363678060&hsa_kw=&hsa_cam=8160495724&gclid=CjwKCAiAgbiQBhAHEiwAuQ6BkpSeIfLXV0Ie_WRTE5Sr3PoflB_ETufbhv945KfX3g-yVb6RIfrBixoC_GIQAvD_BwE>`__ provides an excellent review on the topic.

Mount NTFS Drive on CentOS
==========================

EPEL (Extra Packages for Enterprise Linux) is a Fedora Special Interest Group that creates, maintains, and manages a set of additional high quality packages for Enterprise Linux, including, but not limited to, Red Hat Enterprise Linux (RHEL), CentOS and Scientific Linux (SL), Oracle Linux (OL).

.. code-block:: text

    $ yum install epel-release
    $ yum install ntfs-3g
    
Now your system should be able to open the NTFS drive without issues.

References:

  - `How to Mount a NTFS Drive on CentOS / RHEL / Scientific Linux <https://www.howtoforge.com/tutorial/mount-ntfs-centos/>`__

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

Bomb, Baby!
-----------

`Reference <https://dev.to/itepsilon/foobar-bomb-baby-3l1>`__

`Reference <https://github.com/ivanseed/google-foobar-help/blob/master/challenges/bomb_baby/bomb_baby.md>`__

.. code-block:: note

    def solution(x, y):
        M, F = max(int(x), int(y)), min(int(x), int(y))
        res = 0
        while F > 0:
            res += M // F
            M, F = F, M % F
        if M != 1:
            return 'impossible'
        return str(res - 1)

    print solution('4', '7')         # 4
    print solution('2', '1')         # 1
    print solution('2', '4')         # 'impossible'
    print solution('54000', '30000') # 'impossible'

Fuel Injection Perfection
-------------------------

`Reference <https://gist.github.com/thorstenhirsch/f14842aaeb2d2073e18ec91211ec3875>`__

.. code-block:: note

    def solution(n):
        n = int(n)

        counter = 0

        while n > 3:
            if n & 1:
                if n & 2:
                    n = (n + 1) >> 2
                    counter += 3
                else:
                    n = (n - 1) >> 1
                    counter += 2
            else:
                n = n >> 1
                counter += 1

        if n == 3:
            n = n - 1
            counter += 1

        if n == 2:
            n = n - 1
            counter += 1

        return counter

    print(solution("4"))  # 2
    print(solution("15")) # 5
    
Free the Bunny Workers
----------------------

`Reference <https://vitaminac.github.io/Google-Foobar-Free-the-Bunny-Prisoners/>`__

.. code-block:: note

    from itertools import combinations

    def solution(num_buns, num_required):
        keyrings = [[] for num in range(num_buns)]
        copies_per_key = num_buns - num_required + 1
        for key, bunnies in enumerate(combinations(range(num_buns), copies_per_key)):
            for bunny in bunnies:
                keyrings[bunny].append(key)

        return keyrings

    print solution(2, 1)
    print solution(4, 4)
    print solution(5, 3)
    print solution(3, 1)
    print solution(2, 2)
    print solution(3, 2)
