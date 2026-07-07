/*******************************************************************************
 * Copyright (c) 2017, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.graph.json.test

import com.google.common.base.CharMatcher
import org.eclipse.elk.alg.test.PlainJavaInitialization
import org.eclipse.elk.core.options.CoreOptions
import org.eclipse.elk.core.options.Direction
import org.eclipse.elk.graph.json.ElkGraphJson
import org.eclipse.elk.graph.util.ElkGraphUtil
import org.junit.BeforeClass
import org.junit.Test

import static org.junit.Assert.*

/**
 */
class LayoutOptionsTest {
    
    @BeforeClass
    static def void init() {
        PlainJavaInitialization.initializePlainJavaLayout
    }
    
    @Test
    def void testImportLayoutOptions() {
        val graph = '''
        {
          "id": "root",
          "layoutOptions": {
            "elk.direction": "DOWN"
          }
        }
        '''
        val root = ElkGraphJson.forGraph(graph).toElk
        assertTrue(root.hasProperty(CoreOptions.DIRECTION))
        assertEquals(Direction.DOWN, root.getProperty(CoreOptions.DIRECTION))        
    }

    @Test
    def void testImportPropertiesLegacy() {
        val graph = '''
        {
          "id": "root",
          "properties": {
            "elk.direction": "DOWN"
          }
        }
        '''
        val root = ElkGraphJson.forGraph(graph).toElk
        assertTrue(root.hasProperty(CoreOptions.DIRECTION))
        assertEquals(Direction.DOWN, root.getProperty(CoreOptions.DIRECTION))        
    }    

    @Test
    def void testImportLayoutOptionsHavePriority() {
        val graph = '''
        {
          "id": "root",
          "layoutOptions": {
            "elk.direction": "UP"
          },
          "properties": {
            "elk.direction": "DOWN"
          }
        }
        '''
        val root = ElkGraphJson.forGraph(graph).toElk
        assertTrue(root.hasProperty(CoreOptions.DIRECTION))
        assertEquals(Direction.UP, root.getProperty(CoreOptions.DIRECTION))        
    } 
 
    @Test
    def void exportLayoutOptions() {
        val root = ElkGraphUtil.createGraph
        root.setProperty(CoreOptions.DIRECTION, Direction.UP)
        
        val json = ElkGraphJson.forGraph(root).toJson
        
        val expected = '''{"id":"n0","layoutOptions":{"elk.direction":"UP"}}'''
        val cleanJson = CharMatcher.whitespace.removeFrom(json)
        
        assertEquals(expected, cleanJson)
    }
    
}