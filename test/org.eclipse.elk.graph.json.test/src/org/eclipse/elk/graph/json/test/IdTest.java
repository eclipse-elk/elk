/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
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
package org.eclipse.elk.graph.json.test;

import org.eclipse.elk.graph.json.ElkGraphJson;
import org.eclipse.elk.graph.json.JsonImportException;
import org.junit.Test;

/**
 * Tests for the json format.
 */
public class IdTest {
    
    @Test(expected = JsonImportException.class)
    public void testNoId() {
        ElkGraphJson.forGraph("{}").toElk();
    }
    
    @Test(expected = JsonImportException.class)
    public void testWrongIdTypeNumber() {
        ElkGraphJson.forGraph("{ id: 1.2 }").toElk();
    }

    @Test(expected = JsonImportException.class)
    public void testWrongIdTypeObject() {
        ElkGraphJson.forGraph("{ id: {} }").toElk();
    }

    @Test(expected = JsonImportException.class)
    public void testWrongIdTypeArray() {
        ElkGraphJson.forGraph("{ id: [] }").toElk();
    }
    
    @Test(expected = JsonImportException.class)
    public void testWrongIdTypeBoolean() {
        ElkGraphJson.forGraph("{ id: true }").toElk();
    }
    
    @Test
    public void testGoodIdString() {
        ElkGraphJson.forGraph("{ id: 'foo' }").toElk();
    }
    
    @Test
    public void testGoodIdInt() {
        ElkGraphJson.forGraph("{ id: 3 }").toElk();
    }
    
}
