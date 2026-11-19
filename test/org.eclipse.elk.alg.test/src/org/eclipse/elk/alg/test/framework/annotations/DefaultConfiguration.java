/*******************************************************************************
 * Copyright (c) 2018, 2019 Kiel University and others.
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
package org.eclipse.elk.alg.test.framework.annotations;

import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

/**
 * Instruct the test framework to configure elements of input graphs with default values. The configuration with
 * defaults is done after running configuration methods.
 */
@Retention(RetentionPolicy.RUNTIME)
@Target(ElementType.TYPE)
public @interface DefaultConfiguration {
    
    /** Whether edges should be configured with default values. */
    boolean edges() default true;
    /** Whether nodes should be configured with default values. */
    boolean nodes() default true;
    /** Whether ports should be configured with default values. */
    boolean ports() default true;
    /** Whether random edge labels should be added to the graph. */
    boolean edgeLabels() default false;
    
}
