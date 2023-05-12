# Create SSH keys

```
ssh-keygen -t ed25519
```

It will ask filename and location to save. Please enter following by replacing <akmami> to your local machine username

```
Enter file in which to save the key (/Users/<akmami>/.ssh/id_ed25519): /Users/<akmami>/.ssh/id_ed25519_[server]
```
  
Please do not save the default name since later on, you will forget probably about that key and what server it is being used for etc.

You can enter passphrase if you want or hit enter if not. This helps to prevent others or you (accidentally) to use the public key. In case you use the passphrase, whenever you want to connect to server with ssh, it will ask you the passphrase.

Then, you will see the following image as an output for successful creation of ssh keys.

Now, you have to create config file in order to make your life easier

```
cd ~/.ssh
[# if touch file does not exist]
[touch config]
nano config
```
  
Add following configurations to the config file
 
```
Host [server]
  HostName HOST_NAME
  User USER
  IdentityFile ~/.ssh/id_ed25519_[server].pub
```
  
Also, add new key to your machine, so, it can identify the key 

```
ssh-add -k ~/.ssh/id_ed25519_[server]
```

For some cases or all the time, machine forgot the identity of the keys when the terminal is quit. Hence, each time, you might have to call ssh-add command before, each time you want to connect server.

Now, you have to add ssh key to the server.

You can connect to server with password authentication. After adding the key to the authorized_keys, you can connect to the server by executing following command:

```
ssh [server]
```
