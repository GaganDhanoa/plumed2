#! /bin/bash

actions=$(
eval "$1" manual --action 2>&1 | 
awk '{
  if(NR==1) next;
  if(NF!=1) exit;
  print $1
}' | awk '{
  if($1=="DISTANCE"){
    print $1 ",option:ATOMS,flag:COMPONENTS,flag:SCALED_COMPONENTS,flag:NUMERICAL_DERIVATIVES,flag:NOPBC"
  } else if($1=="MOVINGRESTRAINT"){
    print $1 ",option:ARG,option:VERSE,numbered:STEP,numbered:AT,numbered:KAPPA,flag:NUMERICAL_DERIVATIVES"
  } else print
}'
)

cat << \EOF
" Vim syntax file
" Language: PLUMED

if exists("b:current_syntax")
  finish
endif

" All except space and dash are in word
set iskeyword=33,34,36-126

" Matching label
syntax match   plumedLabel "\v<LABEL\=\S*" contained
highlight link plumedLabel Type

" Matching dots, possibly followed by a comment
" Only dots are part of the match
syntax match plumedDots /\v\.\.\.(\s*(#.*)*$)@=/ contained
highlight link plumedDots Type

EOF
for a in $actions ; do
action_name="${a%%,*}" 
action_name_=$(echo $action_name | sed s/-/_/g)

for l in $(echo "$a" | sed 's/,/ /g')
do
  string=
  case "$l" in
  (flag:*)
# syntax match   plumedKeywordsDISTANCE "\v<COMPONENTS>" contained
    string='"\v<'${l#flag:}'>"' ;;
  (option:*)
# syntax match   plumedKeywordsDISTANCE "\v<ATOMS\=[^{ ]*" contained
    string='"\v<'${l#option:}'\=[^{ #]*"' ;;
  (numbered:*)
# syntax match   plumedKeywordsMOVINGRESTRAINT "\v<KAPPA[0-9]+\=[^{ ]*" contained
    string='"\v<'${l#numbered:}'[0-9]+\=[^{ #]*"' ;;
  esac
  test -n "$string" && echo "syntax match   plumedKeywords$action_name_ $string contained"
done

cat << \EOF | sed s/ACTION/$action_name/g | sed s/ACTNAME/$action_name_/g
" single line, with explicit LABEL
" matching action at beginning of line, till the end of the line
" can contain all the keywords associated with this action, plus strings, label, and comments
syntax region plumedLineACTNAME matchgroup=plumedActionACTNAME start=/\v^\s*ACTION>/ excludenl end=/$/ contains=plumedComment,plumedKeywordsACTNAME,plumedLabel,plumedString fold
" multiple line, with explicit LABEL
" first row might contain extra words before arriving at the dots
" thus continuation dots are matched by plumedDots
" matching action at beginning of line, followed by dots and possibly by a comment
" ends on dots, possibly followed by the same action name and possibly a comment
" comments and initial dots are not part of the match
" can contain all the keywords associated with this action, plus strings, label, and comments
syntax region plumedLineACTNAME matchgroup=plumedActionACTNAME start=/\v^\s*ACTION>(.+\.\.\.\s*(#.*)*$)@=/ end=/\v^\s*\.\.\.(\s+ACTION)?\s*((#.*)*$)@=/ contains=plumedComment,plumedKeywordsACTNAME,plumedLabel,plumedString,plumedDots fold
" single line, with label: syntax
" matching label followed by action
" can contain all the keywords associated with this action, plus strings and comments
" labels are not allwed
syntax region plumedLineACTNAME matchgroup=plumedActionACTNAME start=/\v^\s*\S+:\s+ACTION/ excludenl end=/$/ contains=plumedComment,plumedKeywordsACTNAME,plumedString fold
" multiple line, with label: syntax
" first row might contain extra words before arriving at the dots
" thus continuation dots are matched by plumedDots
" matching label, action, dots, and possibly comment
" comments and dots are not part of the match
" ends on dots, possibly followed by the same label and possibly a comment
" comments and initial dots are not part of the match
syntax region plumedLineACTNAME matchgroup=plumedActionACTNAME start=/\v^\s*\z(\S+\:)\s+ACTION>(.+\.\.\.\s*(#.*)*$)@=/ end=/\v^\s*\.\.\.(\s+\z1)?\s*((#.*)*$)@=/ contains=plumedComment,plumedKeywordsACTNAME,plumedString,plumedDots fold
" this is a hack required to match the ACTION when it is in the second line
syntax match plumedSpecialACTNAME /\v(\.\.\.\s*(#.*)*\_s*)@<=ACTION>/ contained
highlight link plumedSpecialACTNAME Type
" multiple line, with label: syntax
" here ACTION is on the second line
" matching label, dots, possibly comments, newline, then action name
" comments, dots, and action are not part of the match
" ends on dots possibly followed by the same label and possibly a comment
syntax region plumedLineACTNAME matchgroup=plumedActionACTNAME start=/\v^\s*\z(\S+\:)\s+(\.\.\.\s*(#.*)*\_s*ACTION)@=/ end=/\v^\s*\.\.\.(\s+\z1)?\s*((#.*)*$)@=/ contains=plumedComment,plumedKeywordsACTNAME,plumedString,plumedSpecialACTNAME,plumedDots fold
highlight link plumedActionACTNAME Type
highlight link plumedKeywordsACTNAME Statement
EOF

done
cat << \EOF
" comments and strings last, with highest priority
syntax region  plumedString start=/\v\{/  end=/\v\}/
syntax region  plumedString start=/\v\(/  end=/\v\)/
highlight link plumedString String
syntax region  plumedComment start="\v^\S*ENDPLUMED" end="\%$"
syntax match   plumedComment excludenl "\v#.*$"
highlight link plumedComment Comment
EOF

# colors:
# Constant
# Identifier
# Statement
# PreProc
# Type
# Special
# Underlined
# Ignore
# Error
# Todo


