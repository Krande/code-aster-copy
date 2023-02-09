# To contribute to the project, you need to do it through pull/merge request

First you need to fork the repository into your own account. You can do that
simply by clicking the **fork** button on the GitLab interface.

Then, clone the repository on your laptop.

```shell
git clone https://gitlab.com/your-username/your-forkname.git
```

Once this is done, you can setup the *codeaster/src* repository as the upstream
of your clone to simplify the update of your fork repository.

```shell
git remote add upstream https://gitlab.com/codeaster/src.git
```

Now, you have your repository configured, and you want to create a new pull request.
The first step is to create a branch from the HEAD of the **main** branch of
your fork repository.

```shell
git checkout -b your_branch_name
```

Apply your modifications in your branch. Then, you need to push this branch on
your online repository.

```shell
git push -f origin your_branch_name
```

or without `-f`, if the branch already exists online, and you just want to update it.

Once your branch is online, on the GitLab interface, go to the branches webpage,
select the branch you want to push as a merge request, and push the button !

***Be careful to check the 'close after merge' check box, and to push to the
codeaster/src repository***.
