Google Cloud Platform
*********************

`Google Cloud Platform (GCP) <https://cloud.google.com/>`__ is a suite of
cloud computing services that runs on the same infrastructure that Google
uses internally for its end-user products, such as Google Search, Gmail,
Google Drive, and YouTube.

Frequently used commands for GCP
================================

SSH into VM:

.. code-block:: text

    $ gcloud compute ssh <user>@<instance> --project <project> --zone <zone>
    $ gcloud compute ssh "sbslee@instance-1" --project "white-fiber-346900" --zone "us-west1-b"

SSH into VM (username omitted):

.. code-block:: text

    $ gcloud compute ssh <instance> --project <project> --zone <zone>
    $ gcloud compute ssh "instance-1" --project "white-fiber-346900" --zone "us-west1-b"

Upload files to VM:

.. code-block:: text

    $ gcloud compute scp test.txt "instance-1":~ --project "white-fiber-346900" --zone "us-west1-b"
    $ gcloud compute scp --recurse test "instance-1":~ --project "white-fiber-346900" --zone "us-west1-b"

Download files from VM:

.. code-block:: text

    $ gcloud compute scp "instance-1":~/test.txt . --project "white-fiber-346900" --zone "us-west1-b"
    $ gcloud compute scp --recurse "instance-1":~/test . --project "white-fiber-346900" --zone "us-west1-b"

Create VM
=========

- Click ``Create instance``.
- Enter server name.
- Select region (e.g. ``us-central1``).
- Select zone (e.g. ``us-central1-a``)
- Select machine type (e.g. ``e2-micro``).
- Configure boot disk
- Allow HTTP traffic in firewall section.

  * Select OS (e.g. ``Ubuntu``).
  * Select OS version (e.g. ``Ubuntu 18.04 LTS``).
  * Click ``Select``.

- Click ``Create``.
