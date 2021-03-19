# OpenSSH

## Table of Contents

* [Frequently used commands](#Frequently-used-commands)
* [Creating a channel with password](#Creating-a-channel-with-password)
* [Creating a channel without password](#Creating-a-channel-without-password)
* [Channeling through multiple servers](#Channeling-through-multiple-servers)

## Frequently used commands <a name="Frequently-used-commands"></a>

To remove all keys belonging to a host name:

```
ssh-keygen -R host_name
```

To delete a select key from the authentication agent:

```
ssh-add -d ~/.ssh/host_id_rsa.pub
rm ~/.ssh/host_id_rsa
rm ~/.ssh/host_id_rsa.pub
```

## Creating a channel with password <a name="Creating-a-channel-with-password"></a>

First, open your SSH configuration file:

```
vi ~/.ssh/config
```

Next, add the following:

```
Host host_id
    HostName host_name
    User user_name
```

You can now access the server by (you still need to enter your password):

```
ssh host_id
```

## Creating a channel without password <a name="Creating-a-channel-without-password"></a>

First, set up a channel with password as described above. Then, run the following:

```
ssh-keygen -t rsa -b 4096 -C "host_id"
```

Save the private key as `host_id_rsa` and the public key as `host_id_rsa.pub`. Add the private key to the authentication agent:

```
ssh-add ~/.ssh/host_id_rsa
```

Check whether the addition was successful:

```
ssh-add -L
```

Add the public key to the server:

```
cat ~/.ssh/host_id_rsa.pub | ssh host_id 'cat >> ~/.ssh/authorized_keys'
```

Finally, update the configuration:

```
Host host_id
    HostName host_name
    User user_name
    IdentityFile ~/.ssh/host_id_rsa
```

Now, you shouldn't need to enter the password when logging in.

## Channeling through multiple servers <a name="Channeling-through-multiple-servers"></a>

Imagine the server you work on everyday (server C) can only be accessed through another server (server B). Inconveniently, server B can only be accessed through server A. So, your task is to set up a channel that looks like this: local > server A > server B > server C. To do this, you need to set up the SSH configuration as follows:

```
Host host_id_A
    HostName host_name_A
    User user_name_A
    IdentityFile ~/.ssh/host_id_A_rsa

Host host_id_B
    HostName host_name_B
    User user_name_B
    ProxyCommand ssh host_id_A nc %h %p 2> /dev/null
    IdentityFile ~/.ssh/host_id_B_rsa

Host host_id_C
    HostName host_name_C
    User user_name_C
    ProxyCommand ssh host_id_B nc %h %p 2> /dev/null
    IdentityFile ~/.ssh/host_id_C_rsa
```

You can now access server C directly by:

```
ssh host_id_C
```
