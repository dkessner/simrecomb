//
// test_splice.cpp
//
// Copyright 2012 Darren Kessner
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.
//

#include <list>
#include <iostream>
#include <iterator>


using namespace std;


int main()
{
    list<int> l;
    list<int> m;

    for (int i=0; i<10; i++)
    {
        l.push_back(i);
        m.push_back(10+i);
    }

    copy(l.begin(), l.end(), ostream_iterator<int>(cout, " "));
    cout << endl;
    copy(m.begin(), m.end(), ostream_iterator<int>(cout, " "));
    cout << endl;

    const int position = 5;
    list<int>::iterator l_pos = l.begin();
    for (int i=0; i<position; i++) l_pos++;
    list<int>::iterator m_pos = m.begin();
    for (int i=0; i<position; i++) m_pos++;

    cout << "recombining at position " << position << endl;
    list<int> temp;
    temp.splice(temp.begin(), l, l_pos, l.end());
    l.splice(l.end(), m, m_pos, m.end());
    m.splice(m.end(), temp);

    copy(l.begin(), l.end(), ostream_iterator<int>(cout, " "));
    cout << endl;
    copy(m.begin(), m.end(), ostream_iterator<int>(cout, " "));
    cout << endl;
}

