/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p3order.counting;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;

import org.junit.Test;

// CHECKSTYLEOFF MagicNumber
// CHECKSTYLEOFF javadoc
@SuppressWarnings("restriction")
public class BinaryIndexedTreeTest {
    @Test
    public void sumBefore() throws Exception {
        BinaryIndexedTree ft = new BinaryIndexedTree(5);
        ft.add(1);
        ft.add(2);
        ft.add(1);

        assertThat(ft.rank(1), is(0));
        assertThat(ft.rank(2), is(2));
    }

    @Test
    public void size() throws Exception {
        BinaryIndexedTree ft = new BinaryIndexedTree(5);
        ft.add(2);
        ft.add(1);
        ft.add(1);

        assertThat(ft.size(), is(3));
    }

    @Test
    public void removeAll() throws Exception {
        BinaryIndexedTree ft = new BinaryIndexedTree(5);
        ft.add(0);
        ft.add(2);
        ft.add(1);
        ft.add(1);

        ft.removeAll(1);

        assertThat(ft.size(), is(2));
        assertThat(ft.rank(2), is(1));

        ft.removeAll(1);

        assertThat(ft.size(), is(2));
        assertThat(ft.rank(2), is(1));
    }
}
