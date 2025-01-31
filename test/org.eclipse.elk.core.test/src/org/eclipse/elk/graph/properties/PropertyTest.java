/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.graph.properties;

import static org.junit.Assert.assertTrue;

import java.util.EnumSet;
import java.util.List;
import java.util.Set;

import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.math.ElkMargin;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.options.SizeConstraint;
import org.junit.BeforeClass;
import org.junit.Test;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

/**
 * Tests a set of commonly used {@link Property}s, e.g. that they return a proper default value. This test is not part
 * of the {@code org.eclipse.elk.graph.test} plug-in due to the fact that it needs core classes to work.
 */
public class PropertyTest {

    @BeforeClass
    public static void initMetaDataService() {
        // init
        LayoutMetaDataService.getInstance();
    }
    
    @Test
    public void testPropertyDefaultPrimitive() {
        int i = 43;
        testPropertyPrimitive(i);
        
        float f = 23.3f;
        testPropertyPrimitive(f);
        
        double d = 32.3d;
        testPropertyPrimitive(d);

        String s = "foo";
        testPropertyPrimitive(s);
    }
        
    @Test
    public void testPropertyDefaultIDataObject() {
        // all IDataObjects 
        KVector v = new KVector(2, 3);
        testPropertyObject(v);
        
        KVectorChain vc = new KVectorChain(v, v);
        testPropertyObject(vc);
        
        ElkPadding ep = new ElkPadding(2, 3);
        testPropertyObject(ep);
        
        ElkMargin em = new ElkMargin(3, 2);
        testPropertyObject(em);
    }
    
    @Test
    public void testPropertyDefaultObject() {
        KVector v = new KVector(3, 2);
        
        List<KVector> al = Lists.newArrayList(v, v.normalize(), v.negate());
        testPropertyObject(al);
        
        List<KVector> ll = Lists.newLinkedList(al);
        testPropertyObject(ll);
        
        Set<KVector> hs = Sets.newHashSet(v, v.normalize(), v.negate());
        testPropertyObject(hs);

        Set<KVector> lhs = Sets.newLinkedHashSet(al);
        testPropertyObject(lhs);
        
        Set<Integer> ts = Sets.newTreeSet(Lists.newArrayList(1, 2, 3, 4));
        testPropertyObject(ts);
    }

    @Test
    public void testPropertyDefaultEnum() {
        PortSide ps = PortSide.EAST;
        testPropertyPrimitive(ps);
    }
    
    @Test
    public void testPropertyDefaultEnumSet() {
        EnumSet<SizeConstraint> sc = EnumSet.allOf(SizeConstraint.class);
        testPropertyObject(sc);
    }
    
    private <T> void testPropertyObject(T defaultValue) {
        IProperty<T> p = new Property<T>("dummy" + defaultValue.getClass().getSimpleName(), defaultValue);
        T copy = p.getDefault();
        assertTrue("Default value is not a copy: " + defaultValue.getClass(), defaultValue != copy);
        assertTrue("Copied value is not equal: " + defaultValue.getClass(), defaultValue.equals(copy));
    }
    
    private <T> void testPropertyPrimitive(T defaultValue) {
        IProperty<T> p = new Property<T>("dummy" + defaultValue.getClass().getSimpleName(), defaultValue);
        T copy = p.getDefault();
        assertTrue("Copied value is not equal: " + defaultValue.getClass(), defaultValue == copy);
        assertTrue("Copied value is not equal: " + defaultValue.getClass(), defaultValue.equals(copy));
    }
    
    @Test(expected = IllegalStateException.class)
    public void testUnknownPropertyGetDefault() {
        Cloneable o = new Cloneable() { };
        IProperty<Object> p = new Property<>("o", o);
        p.getDefault();
    }
    
}
