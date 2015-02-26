### reset repo to state at a certain commit:

git reset --hard <commit-hash>
git push -f origin master

### merge changes to a single file rather than commits:

Need to merge just file f of branch B into file f of branch A
(assumes that all changes are committed in both branches A and B):

git checkout A
git checkout --patch B f

see http://stackoverflow.com/a/11593308/295025
