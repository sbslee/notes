Google Cloud Platform
*********************

`Google Cloud Platform (GCP) <https://cloud.google.com/>`__ is a suite of
cloud computing services that runs on the same infrastructure that Google
uses internally for its end-user products, such as Google Search, Gmail,
Google Drive, and YouTube.

Frequently used commands for GCP
================================

SSH into an instance:

.. code-block:: text

    $ gcloud compute ssh --zone "us-west1-b" "instance-1" --project "white-fiber-346900"

Upload files to an instance:

.. code-block:: text

    $ gcloud compute scp test.txt "instance-1":~ --project "white-fiber-346900" --zone "us-west1-b"

Download files from an instance:

.. code-block:: text

    $ gcloud compute scp "instance-1":~/test.txt . --project "white-fiber-346900" --zone "us-west1-b"
